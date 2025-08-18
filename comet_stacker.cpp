#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>

// Qt headers
#include <QCoreApplication>
#include <QString>
#include <QDir>
#include <QFileInfo>
#include <QDateTime>

// OpenCV for image processing and transformations
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>

// CFITSIO for reading FITS files
#include <fitsio.h>

const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;

struct CometPosition {
    double jd;
    double ra;   // degrees
    double dec;  // degrees
};

struct ImageData {
    QString filename;
    double jd;
    double ra_center;    // From WCS
    double dec_center;   // From WCS
    double pixel_scale;  // arcsec/pixel
    double rotation;     // Image rotation (degrees)
    cv::Mat image;
    cv::Size size;
    
    // Calculated comet position
    double comet_ra;
    double comet_dec;
    double comet_x_pixel;
    double comet_y_pixel;
};

class RotationCorrectedCometStacker {
private:
    std::vector<CometPosition> ephemeris;
    std::vector<ImageData> images;
    QString ephemeris_file;
    
public:
    bool loadEphemeris(const QString& filename) {
        ephemeris_file = filename;
        std::ifstream file(filename.toStdString());
        if (!file.is_open()) {
            std::cout << "Cannot open ephemeris file: " << filename.toStdString() << std::endl;
            return false;
        }
        
        std::cout << "Loading ephemeris from: " << filename.toStdString() << std::endl;
        
        std::string line;
        int count = 0;
        
        // Skip header lines
        while (std::getline(file, line) && line[0] == '#') {}
        
        do {
            if (line.empty()) continue;
            
            std::istringstream ss(line);
            std::string jd_str, ra_str, dec_str;
            
            if (std::getline(ss, jd_str, ',') && 
                std::getline(ss, ra_str, ',') && 
                std::getline(ss, dec_str, ',')) {
                
                try {
                    CometPosition pos;
                    pos.jd = std::stod(jd_str);
                    pos.ra = std::stod(ra_str);
                    pos.dec = std::stod(dec_str);
                    ephemeris.push_back(pos);
                    count++;
                } catch (...) {
                    // Skip invalid lines
                }
            }
        } while (std::getline(file, line));
        
        std::cout << "Loaded " << count << " ephemeris points" << std::endl;
        return count > 0;
    }
    
    CometPosition interpolateComet(double jd) {
        if (ephemeris.empty()) {
            return {jd, 0.0, 0.0};
        }
        
        // Find bracketing points
        auto upper = std::lower_bound(ephemeris.begin(), ephemeris.end(), jd,
            [](const CometPosition& pos, double jd) { return pos.jd < jd; });
        
        if (upper == ephemeris.begin()) {
            return ephemeris.front();
        }
        if (upper == ephemeris.end()) {
            return ephemeris.back();
        }
        
        auto lower = upper - 1;
        
        // Linear interpolation
        double t = (jd - lower->jd) / (upper->jd - lower->jd);
        
        CometPosition result;
        result.jd = jd;
        result.ra = lower->ra + t * (upper->ra - lower->ra);
        result.dec = lower->dec + t * (upper->dec - lower->dec);
        
        // Handle RA wraparound
        double ra_diff = upper->ra - lower->ra;
        if (ra_diff > 180.0) {
            result.ra = lower->ra + t * (ra_diff - 360.0);
        } else if (ra_diff < -180.0) {
            result.ra = lower->ra + t * (ra_diff + 360.0);
        }
        
        if (result.ra < 0) result.ra += 360.0;
        if (result.ra >= 360.0) result.ra -= 360.0;
        
        return result;
    }
    
    bool loadImages(const QStringList& filenames) {
        std::cout << "\n=== LOADING IMAGES ===" << std::endl;
        
        for (const QString& filename : filenames) {
            ImageData img;
            if (loadFITSImage(filename, img)) {
                
                // Get comet position for this image time
                auto comet_pos = interpolateComet(img.jd);
                img.comet_ra = comet_pos.ra;
                img.comet_dec = comet_pos.dec;
                
                // Convert comet RA/Dec to pixel coordinates
                convertRaDecToPixel(img, comet_pos.ra, comet_pos.dec, 
                                   img.comet_x_pixel, img.comet_y_pixel);
                
                std::cout << "Image " << images.size() + 1 << ": " 
                          << QFileInfo(filename).fileName().toStdString() << std::endl;
                std::cout << "  JD: " << img.jd << std::endl;
                std::cout << "  Comet predicted: RA=" << comet_pos.ra << "°, Dec=" << comet_pos.dec << "°" << std::endl;
                std::cout << "  Comet pixel: (" << img.comet_x_pixel << ", " << img.comet_y_pixel << ")" << std::endl;
                
                images.push_back(img);
            }
        }
        
        std::cout << "Successfully loaded " << images.size() << " images" << std::endl;
        return !images.empty();
    }
    
private:
    bool loadFITSImage(const QString& filename, ImageData& img) {
        img.filename = filename;
        
        fitsfile *fptr = nullptr;
        int status = 0;
        
        QByteArray pathBytes = filename.toLocal8Bit();
        if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
            std::cout << "Failed to open: " << filename.toStdString() << std::endl;
            return false;
        }
        
        // Read WCS and time information
        double crval1, crval2, cdelt1, cdelt2, crota1 = 0.0;
        char dateobs[FLEN_VALUE];
        
        fits_read_key(fptr, TDOUBLE, "CRVAL1", &crval1, nullptr, &status);
        fits_read_key(fptr, TDOUBLE, "CRVAL2", &crval2, nullptr, &status);
        fits_read_key(fptr, TDOUBLE, "CDELT1", &cdelt1, nullptr, &status);
        fits_read_key(fptr, TDOUBLE, "CDELT2", &cdelt2, nullptr, &status);
        fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, nullptr, &status);
        
        // Try to read rotation information
        int rot_status = 0;
        fits_read_key(fptr, TDOUBLE, "CROTA1", &crota1, nullptr, &rot_status);
        if (rot_status != 0) {
            fits_read_key(fptr, TDOUBLE, "CROTA2", &crota1, nullptr, &rot_status);
        }
        
        if (status != 0) {
            fits_close_file(fptr, &status);
            std::cout << "Failed to read WCS from: " << filename.toStdString() << std::endl;
            return false;
        }
        
        img.ra_center = crval1;
        img.dec_center = crval2;
        img.pixel_scale = abs(cdelt1) * 3600.0; // Convert to arcsec/pixel
        img.rotation = crota1; // Image rotation in degrees
        
        // Parse observation time
        QString dateStr = QString::fromLatin1(dateobs).trimmed().remove('\'').remove('"');
        QDateTime obsTime = QDateTime::fromString(dateStr, "yyyy-MM-ddThh:mm:ss");
        if (!obsTime.isValid()) {
            obsTime = QDateTime::fromString(dateStr, "yyyy-MM-ddThh:mm:ss.zzz");
        }
        
        // Convert to Julian Date (simplified)
        QDate date = obsTime.date();
        QTime time = obsTime.time();
        double fracDay = (time.hour() + time.minute()/60.0 + time.second()/3600.0) / 24.0;
        img.jd = gregorianToJD(date.year(), date.month(), date.day()) + fracDay;
        
        // Read image data
        long naxes[2];
        fits_get_img_size(fptr, 2, naxes, &status);
        img.size = cv::Size(naxes[0], naxes[1]);
        
        long totalPixels = naxes[0] * naxes[1];
        std::vector<float> pixels(totalPixels);
        fits_read_img(fptr, TFLOAT, 1, totalPixels, nullptr, pixels.data(), nullptr, &status);
        
        img.image = cv::Mat(naxes[1], naxes[0], CV_32F, pixels.data()).clone();
        
        fits_close_file(fptr, &status);
        
        return status == 0;
    }
    
    double gregorianToJD(int year, int month, int day) {
        if (month <= 2) {
            year -= 1;
            month += 12;
        }
        int a = year / 100;
        int b = 2 - a + (a / 4);
        return floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + b - 1524.5;
    }
    
    void convertRaDecToPixel(const ImageData& img, double ra, double dec, double& x, double& y) {
        // Simple linear WCS transformation (could be improved with full SIP/distortion)
        double delta_ra = (ra - img.ra_center) * cos(dec * DEG2RAD);
        double delta_dec = dec - img.dec_center;
        
        // Convert to pixels (assuming CRPIX is image center)
        x = img.size.width / 2.0 + delta_ra / (img.pixel_scale / 3600.0);
        y = img.size.height / 2.0 + delta_dec / (img.pixel_scale / 3600.0);
    }
    
public:
    bool stackImages(const QString& output_filename) {
        if (images.empty()) {
            std::cout << "No images to stack" << std::endl;
            return false;
        }
        
        std::cout << "\n=== ROTATION-CORRECTED COMET STACKING ===" << std::endl;
        
        // Use first image as reference
        const ImageData& reference = images[0];
        cv::Mat stacked = cv::Mat::zeros(reference.size, CV_32F);
        cv::Mat weight_map = cv::Mat::zeros(reference.size, CV_32F);
        
        std::cout << "Reference image: " << QFileInfo(reference.filename).fileName().toStdString() << std::endl;
        std::cout << "Reference comet position: (" << reference.comet_x_pixel 
                  << ", " << reference.comet_y_pixel << ")" << std::endl;
        std::cout << "Reference rotation: " << reference.rotation << "°" << std::endl;
        
        for (size_t i = 0; i < images.size(); ++i) {
            const ImageData& img = images[i];
            
            // Calculate transformation to align comet positions with rotation correction
            cv::Point2f ref_center(reference.comet_x_pixel, reference.comet_y_pixel);
            cv::Point2f img_center(img.comet_x_pixel, img.comet_y_pixel);
            
            // Calculate rotation difference
            double rotation_diff = reference.rotation - img.rotation;
            
            // Create transformation matrix
            cv::Mat transform = calculateTransform(img_center, ref_center, rotation_diff, img.size);
            
            std::cout << "Image " << (i+1) << " transformation:" << std::endl;
            std::cout << "  Translation: (" << (ref_center.x - img_center.x) 
                      << ", " << (ref_center.y - img_center.y) << ")" << std::endl;
            std::cout << "  Rotation: " << rotation_diff << "°" << std::endl;
            
            // Apply transformation and accumulate
            cv::Mat transformed;
            cv::warpAffine(img.image, transformed, transform, reference.size, 
                          cv::INTER_LINEAR, cv::BORDER_CONSTANT, cv::Scalar(0));
            
            // Accumulate with proper weighting
            cv::Mat mask = (transformed > 0);
            mask.convertTo(mask, CV_32F, 1.0/255.0);
            
            stacked += transformed;
            weight_map += mask;
        }
        
        // Normalize by weight
        cv::Mat normalized;
        cv::divide(stacked, weight_map, normalized, 1.0, CV_32F);
        
        // Handle areas with no coverage
        normalized.setTo(0, weight_map == 0);
        
        return saveFITSImage(normalized, output_filename, reference);
    }
    
private:
    cv::Mat calculateTransform(cv::Point2f src_center, cv::Point2f dst_center, 
                              double rotation_degrees, cv::Size img_size) {
        
        // Image center for rotation
        cv::Point2f image_center(img_size.width / 2.0f, img_size.height / 2.0f);
        
        // Create rotation matrix around image center
        cv::Mat rotation_matrix = cv::getRotationMatrix2D(image_center, rotation_degrees, 1.0);
        
        // Add translation to align comet centers
        cv::Point2f translation = dst_center - src_center;
        rotation_matrix.at<double>(0, 2) += translation.x;
        rotation_matrix.at<double>(1, 2) += translation.y;
        
        return rotation_matrix;
    }
    
    bool saveFITSImage(const cv::Mat& image, const QString& filename, const ImageData& reference) {
        fitsfile *fptr = nullptr;
        int status = 0;
        
        QByteArray pathBytes = QString("!%1").arg(filename).toLocal8Bit();
        if (fits_create_file(&fptr, pathBytes.data(), &status)) {
            return false;
        }
        
        long naxes[2] = {image.cols, image.rows};
        fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);
        
        // Convert image data
        std::vector<float> pixels(image.rows * image.cols);
        for (int y = 0; y < image.rows; ++y) {
            for (int x = 0; x < image.cols; ++x) {
                pixels[y * image.cols + x] = image.at<float>(y, x);
            }
        }
        
        fits_write_img(fptr, TFLOAT, 1, pixels.size(), pixels.data(), &status);
        
        // Copy WCS from reference image
        fits_write_key(fptr, TDOUBLE, "CRVAL1", (void *)&reference.ra_center, "Reference RA", &status);
        fits_write_key(fptr, TDOUBLE, "CRVAL2", (void *)&reference.dec_center, "Reference Dec", &status);
        
        double crpix1 = reference.size.width / 2.0;
        double crpix2 = reference.size.height / 2.0;
        fits_write_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "Reference pixel X", &status);
        fits_write_key(fptr, TDOUBLE, "CRPIX2", &crpix2, "Reference pixel Y", &status);
        
        double cdelt1 = -reference.pixel_scale / 3600.0;
        double cdelt2 = reference.pixel_scale / 3600.0;
        fits_write_key(fptr, TDOUBLE, "CDELT1", &cdelt1, "Pixel scale X", &status);
        fits_write_key(fptr, TDOUBLE, "CDELT2", &cdelt2, "Pixel scale Y", &status);
        
        fits_write_key(fptr, TDOUBLE, "CROTA1", (void *)&reference.rotation, "Image rotation", &status);
        
        char ctype1[] = "RA---TAN";
        char ctype2[] = "DEC--TAN";
        fits_write_key(fptr, TSTRING, "CTYPE1", ctype1, "Coordinate type", &status);
        fits_write_key(fptr, TSTRING, "CTYPE2", ctype2, "Coordinate type", &status);
        
        // Add stacking info
        int nimages = images.size();
        fits_write_key(fptr, TINT, "NSTACKED", &nimages, "Number of stacked images", &status);
        
        char method[] = "COMET_ROTATION_CORRECTED";
        fits_write_key(fptr, TSTRING, "STACKMET", method, "Stacking method", &status);
        
        fits_close_file(fptr, &status);
        
        std::cout << "Rotation-corrected stacked image saved: " << filename.toStdString() << std::endl;
        return status == 0;
    }
};

int main(int argc, char* argv[]) {
    QCoreApplication app(argc, argv);
    
    if (argc < 4) {
        std::cout << "Rotation-Corrected Comet Stacker" << std::endl;
        std::cout << "Usage: " << argv[0] << " <ephemeris.csv> <image_dir> <output.fits> [max_images]" << std::endl;
        std::cout << std::endl;
        std::cout << "This version corrects for field rotation during stacking." << std::endl;
        std::cout << "Example: " << argv[0] << " ephemeris.csv /path/to/lights comet_stack_corrected.fits 100" << std::endl;
        return 1;
    }
    
    QString ephemeris_file = QString::fromStdString(argv[1]);
    QString image_dir = QString::fromStdString(argv[2]);
    QString output_file = QString::fromStdString(argv[3]);
    int max_images = (argc > 4) ? std::stoi(argv[4]) : 100;
    
    RotationCorrectedCometStacker stacker;
    
    // Load ephemeris
    if (!stacker.loadEphemeris(ephemeris_file)) {
        std::cout << "Failed to load ephemeris file" << std::endl;
        return 1;
    }
    
    // Find FITS files
    QDir dir(image_dir);
    QStringList filters;
    filters << "*.fits" << "*.fit" << "*.FITS" << "*.FIT";
    QStringList fits_files = dir.entryList(filters, QDir::Files, QDir::Name);
    
    if (fits_files.isEmpty()) {
        std::cout << "No FITS files found in: " << image_dir.toStdString() << std::endl;
        return 1;
    }
    
    // Limit to first N images
    if (fits_files.size() > max_images) {
        fits_files = fits_files.mid(0, max_images);
    }
    
    // Convert to full paths
    QStringList full_paths;
    for (const QString& file : fits_files) {
        full_paths.append(dir.absoluteFilePath(file));
    }
    
    std::cout << "Found " << fits_files.size() << " FITS files (limiting to " << max_images << ")" << std::endl;
    
    // Load and stack images
    if (stacker.loadImages(full_paths)) {
        if (stacker.stackImages(output_file)) {
            std::cout << "\n✓ Rotation-corrected comet stack completed successfully!" << std::endl;
            std::cout << "Output: " << output_file.toStdString() << std::endl;
            std::cout << "\nStars should now appear as points instead of curved trails." << std::endl;
        }
    }
    
    return 0;
}
