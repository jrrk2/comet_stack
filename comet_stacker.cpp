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
#include <QStringList>
#include <QJsonDocument>
#include <QJsonObject>
#include <QFile>
#include <QRegularExpression>

// OpenCV for image processing and transformations
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>

// CFITSIO for reading FITS files
#include <fitsio.h>

// PCL for robust WCS handling
#include "StarCatalogValidator.h"
#include "ImageReader.h"

const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;
const double SIDEREAL_RATE = 15.041067179; // arcsec/second

enum class TrackingMode {
    SIDEREAL,    // Track at sidereal rate (stars stationary)
    COMET,       // Track comet motion (comet stationary)
    FIXED        // No tracking (everything moves)
};

struct StackingConfig {
    TrackingMode mode = TrackingMode::SIDEREAL;
    
    // Sidereal tracking parameters
    double target_ra = 0.0;    // Target RA for sidereal tracking (degrees)
    double target_dec = 0.0;   // Target Dec for sidereal tracking (degrees)
    double reference_jd = 0.0; // Reference time for tracking
    
    // Comet tracking parameters
    QString ephemeris_file;
    
    // Image processing
    int max_images = 100;
    double rotation_tolerance = 0.1; // degrees
    bool apply_rotation_correction = true;
    
    // Output settings
    QString output_file = "stacked.fits";
    bool save_diagnostics = false;
    
    void printConfig() {
        std::cout << "\n=== STACKING CONFIGURATION ===" << std::endl;
        
        switch(mode) {
            case TrackingMode::SIDEREAL:
                std::cout << "Mode: SIDEREAL TRACKING" << std::endl;
                std::cout << "Target: RA=" << target_ra << "° Dec=" << target_dec << "°" << std::endl;
                std::cout << "Reference JD: " << reference_jd << std::endl;
                break;
            case TrackingMode::COMET:
                std::cout << "Mode: COMET TRACKING" << std::endl;
                std::cout << "Ephemeris: " << ephemeris_file.toStdString() << std::endl;
                break;
            case TrackingMode::FIXED:
                std::cout << "Mode: FIXED (no tracking)" << std::endl;
                break;
        }
        
        std::cout << "Max images: " << max_images << std::endl;
        std::cout << "Rotation correction: " << (apply_rotation_correction ? "ON" : "OFF") << std::endl;
        std::cout << "Output: " << output_file.toStdString() << std::endl;
    }
};

struct CometPosition {
    double jd;
    double ra;   // degrees
    double dec;  // degrees
};

struct EnhancedImageData {
    QString filename;
    double jd;
    cv::Mat image;
    cv::Mat rawImage;  // Original image before dark subtraction
    cv::Size size;
    
    // Enhanced WCS data (using PCL)
    std::unique_ptr<StarCatalogValidator> wcsValidator;
    bool hasValidWCS;
    
    // Target tracking positions
    QPointF target_pixel;    // Where we want to track (sidereal point or comet)
    double pixel_scale;      // arcsec/pixel  
    QPointF image_center_radec;
    
    // Image metadata for dark matching
    int exposure_time = 10;     // seconds
    int temperature = 273;      // Kelvin
    QString binning = "1x1";
    QString bayer_pattern = "RGGB";
    
    // Dark frame processing
    bool darkSubtracted = false;
    QString matchedDarkFrame;
    
    // Diagnostics
    double calculated_rotation = 0.0;
    QPointF translation_offset;
    
    EnhancedImageData() : hasValidWCS(false), pixel_scale(0.0), darkSubtracted(false) {}
};

class ConfigurableStacker {
private:
    StackingConfig config;
    std::vector<CometPosition> ephemeris;
    std::vector<EnhancedImageData> images;
    QString m_darkDirectory;
    
    // Dark frame processing parameters
    double m_exposureTolerance = 25.0;  // 25% tolerance (more permissive)
    int m_temperatureTolerance = 15;    // 15K tolerance (more permissive)
    
public:
    void setConfig(TrackingMode mode, double ra = 0.0, double dec = 0.0, QString ephemeris = "", QString darkDir = "") {
        config.mode = mode;
        config.target_ra = ra;
        config.target_dec = dec;
        config.ephemeris_file = ephemeris;
        m_darkDirectory = darkDir;
    }
    
    bool loadEphemeris() {
        if (config.mode != TrackingMode::COMET || config.ephemeris_file.isEmpty()) {
            return true; // Not needed for sidereal mode
        }
        
        std::ifstream file(config.ephemeris_file.toStdString());
        if (!file.is_open()) {
            std::cout << "Cannot open ephemeris file: " << config.ephemeris_file.toStdString() << std::endl;
            return false;
        }
        
        std::cout << "Loading ephemeris from: " << config.ephemeris_file.toStdString() << std::endl;
        
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
        
        std::cout << "Loaded " << count << " observations" << std::endl;
        return count > 0;
    }
    
    QPointF calculateSiderealTarget(double jd, const EnhancedImageData& img) {
        // Calculate where the target RA/Dec should be at this time
        // accounting for sidereal motion
        
        double time_diff_hours = (jd - config.reference_jd) * 24.0;
        double sidereal_motion_arcsec = time_diff_hours * SIDEREAL_RATE;
        double sidereal_motion_deg = sidereal_motion_arcsec / 3600.0;
        
        // The target moves westward due to Earth's rotation
        double target_ra_now = config.target_ra - sidereal_motion_deg;
        
        // Normalize RA
        while (target_ra_now < 0) target_ra_now += 360.0;
        while (target_ra_now >= 360.0) target_ra_now -= 360.0;
        
        // Convert to pixel coordinates
        return img.wcsValidator->skyToPixel(target_ra_now, config.target_dec);
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
        config.printConfig();
        
        int failed_wcs_count = 0;
        
        for (const QString& filename : filenames) {
            EnhancedImageData img;
            
            if (loadFITSImageWithPCL(filename, img)) {
                
                // Calculate target position based on tracking mode
                switch(config.mode) {
                    case TrackingMode::SIDEREAL:
                        img.target_pixel = calculateSiderealTarget(img.jd, img);
                        std::cout << "Image " << images.size() + 1 << ": " 
                                  << QFileInfo(filename).fileName().toStdString() << std::endl;
                        std::cout << "  JD: " << img.jd << std::endl;
                        std::cout << "  Sidereal target pixel: (" << img.target_pixel.x() 
                                  << ", " << img.target_pixel.y() << ")" << std::endl;
                        break;
                        
                    case TrackingMode::COMET: {
                        auto comet_pos = interpolateComet(img.jd);
                        img.target_pixel = img.wcsValidator->skyToPixel(comet_pos.ra, comet_pos.dec);
                        std::cout << "Image " << images.size() + 1 << ": " 
                                  << QFileInfo(filename).fileName().toStdString() << std::endl;
                        std::cout << "  JD: " << img.jd << std::endl;
                        std::cout << "  Comet: RA=" << comet_pos.ra << "° Dec=" << comet_pos.dec << "°" << std::endl;
                        std::cout << "  Comet pixel: (" << img.target_pixel.x() 
                                  << ", " << img.target_pixel.y() << ")" << std::endl;
                        break;
                    }
                    
                    case TrackingMode::FIXED:
                        img.target_pixel = QPointF(img.size.width/2.0, img.size.height/2.0);
                        std::cout << "Image " << images.size() + 1 << ": " 
                                  << QFileInfo(filename).fileName().toStdString() << std::endl;
                        std::cout << "  JD: " << img.jd << std::endl;
                        std::cout << "  Fixed center: (" << img.target_pixel.x() 
                                  << ", " << img.target_pixel.y() << ")" << std::endl;
                        break;
                }
                
                images.push_back(std::move(img));
            } else {
                failed_wcs_count++;
            }
        }
        
        std::cout << "Successfully loaded " << images.size() << " images" << std::endl;
        if (failed_wcs_count > 0) {
            std::cout << "Failed to load WCS from " << failed_wcs_count << " images" << std::endl;
        }
        
        return !images.empty();
    }
    
private:
    bool loadFITSImageWithPCL(const QString& filename, EnhancedImageData& img) {
        img.filename = filename;
        img.wcsValidator = std::make_unique<StarCatalogValidator>();
        
        fitsfile *fptr = nullptr;
        int status = 0;
        
        QByteArray pathBytes = filename.toLocal8Bit();
        if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
            std::cout << "Failed to open: " << filename.toStdString() << std::endl;
            return false;
        }
        
        // Read basic image properties
        long naxes[2];
        fits_get_img_size(fptr, 2, naxes, &status);
        img.size = cv::Size(naxes[0], naxes[1]);
        
        // Read observation time
        char dateobs[FLEN_VALUE];
        status = 0;
        if (fits_read_key(fptr, TSTRING, "DATE-OBS", dateobs, nullptr, &status) == 0) {
            QString dateStr = QString::fromLatin1(dateobs).trimmed().remove('\'').remove('"');
            QDateTime obsTime = QDateTime::fromString(dateStr, "yyyy-MM-ddThh:mm:ss");
            if (!obsTime.isValid()) {
                obsTime = QDateTime::fromString(dateStr, "yyyy-MM-ddThh:mm:ss.zzz");
            }
            
            if (obsTime.isValid()) {
                // Convert to Julian Date
                QDate date = obsTime.date();
                QTime time = obsTime.time();
                double fracDay = (time.hour() + time.minute()/60.0 + time.second()/3600.0) / 24.0;
                img.jd = gregorianToJD(date.year(), date.month(), date.day()) + fracDay;
            }
        }
        
        // Create ImageData structure for PCL WCS setup
        ImageData imageData;
        imageData.width = naxes[0];
        imageData.height = naxes[1];
        
        // Read all FITS keywords for PCL
        int nkeys;
        fits_get_hdrspace(fptr, &nkeys, nullptr, &status);
        
        for (int i = 1; i <= nkeys; ++i) {
            char card[FLEN_CARD];
            if (fits_read_record(fptr, i, card, &status) == 0) {
                QString cardStr = QString::fromLatin1(card);
                
                // Parse FITS card into key-value format that PCL expects
                QRegularExpression fitsRegex(R"(^([A-Z0-9_-]{1,8})\s*=\s*([^/]*?)(?:\s*/\s*(.*))?$)");
                auto match = fitsRegex.match(cardStr);
                
                if (match.hasMatch()) {
                    QString key = match.captured(1).trimmed();
                    QString value = match.captured(2).trimmed();
                    QString comment = match.captured(3).trimmed();
                    
                    // Remove quotes from string values
                    if (value.startsWith("'") && value.endsWith("'")) {
                        value = value.mid(1, value.length() - 2).trimmed();
                    }
                    
                    QString metadataLine = QString("%1: %2").arg(key).arg(value);
                    if (!comment.isEmpty()) {
                        metadataLine += QString(" (%1)").arg(comment);
                    }
                    
                    imageData.metadata.append(metadataLine);
                }
            }
        }
        
        // Setup PCL WCS using enhanced validator
        img.hasValidWCS = img.wcsValidator->setWCSFromImageMetadata(imageData);
        
        if (img.hasValidWCS) {
            // Extract diagnostic information using PCL
            img.pixel_scale = img.wcsValidator->getPixScale();
            
            double centerRA, centerDec;
            if (img.wcsValidator->getCenter(centerRA, centerDec)) {
                img.image_center_radec = QPointF(centerRA, centerDec);
            }
        }
        
        // Read image data
        long totalPixels = naxes[0] * naxes[1];
        std::vector<float> pixels(totalPixels);
        fits_read_img(fptr, TFLOAT, 1, totalPixels, nullptr, pixels.data(), nullptr, &status);
        
        // Store raw image data
        img.rawImage = cv::Mat(naxes[1], naxes[0], CV_32F, pixels.data()).clone();
        
        // Extract metadata for dark frame matching
        extractImageMetadata(fptr, img);
        
        fits_close_file(fptr, &status);
        
        if (status != 0) {
            return false;
        }
        
        // Apply dark frame subtraction if dark directory is specified
        if (!m_darkDirectory.isEmpty()) {
            applyDarkFrameSubtraction(img);
        } else {
            // If no dark frames, use raw image
            img.image = img.rawImage.clone();
        }
        
        return img.hasValidWCS;
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
    
    // Dark frame processing methods (adapted from StellinaProcessor)
    void extractImageMetadata(fitsfile* fptr, EnhancedImageData& img) {
        int status = 0;
        
        // Extract exposure time (try multiple keywords and units)
        double exposure = 10.0;
        if (fits_read_key(fptr, TDOUBLE, "EXPTIME", &exposure, nullptr, &status) == 0) {
            // EXPTIME usually in seconds
            img.exposure_time = static_cast<int>(exposure);
        } else if (fits_read_key(fptr, TDOUBLE, "EXPOSURE", &exposure, nullptr, &status) == 0) {
            // EXPOSURE might be in milliseconds for Stellina
            if (exposure > 100) {
                // Likely milliseconds, convert to seconds
                img.exposure_time = static_cast<int>(exposure / 1000.0);
                std::cout << "  Exposure: " << exposure << "ms -> " << img.exposure_time << "s" << std::endl;
            } else {
                // Already in seconds
                img.exposure_time = static_cast<int>(exposure);
            }
        }
        status = 0;
        
        // Extract temperature (handle Celsius to Kelvin conversion)
        double temp = 15.0; // Default 15°C
        if (fits_read_key(fptr, TDOUBLE, "TEMP", &temp, nullptr, &status) == 0 ||
            fits_read_key(fptr, TDOUBLE, "CCD-TEMP", &temp, nullptr, &status) == 0 ||
            fits_read_key(fptr, TDOUBLE, "TEMPERAT", &temp, nullptr, &status) == 0) {
            
            // For Stellina, TEMP is in Celsius, convert to Kelvin
            if (temp < 100) {
                temp += 273.15;
                std::cout << "  Temperature: " << (temp - 273.15) << "°C -> " << temp << "K" << std::endl;
            }
            img.temperature = static_cast<int>(temp);
        }
        status = 0;
        
        // Extract binning
        char binning[FLEN_VALUE] = "1x1";
        if (fits_read_key(fptr, TSTRING, "XBINNING", binning, nullptr, &status) == 0) {
            img.binning = QString::fromLatin1(binning).trimmed().remove('\'').remove('"');
        }
        status = 0;
        
        // Extract bayer pattern
        char bayer[FLEN_VALUE] = "RGGB";
        if (fits_read_key(fptr, TSTRING, "BAYERPAT", bayer, nullptr, &status) == 0 ||
            fits_read_key(fptr, TSTRING, "COLORTYP", bayer, nullptr, &status) == 0) {
            img.bayer_pattern = QString::fromLatin1(bayer).trimmed().remove('\'').remove('"').toUpper();
        }
        
        // Check if this is a reversed Stellina image (img-####r pattern)
        QString basename = QFileInfo(img.filename).baseName();
        if (basename.contains(QRegularExpression(R"(img-\d+r$)", QRegularExpression::CaseInsensitiveOption))) {
            std::cout << "  Detected reversed Stellina image - using BGGR pattern" << std::endl;
            img.bayer_pattern = "BGGR";
        }
        
        std::cout << "  Metadata: " << img.exposure_time << "s, " << img.temperature << "K, " 
                  << img.binning.toStdString() << ", " << img.bayer_pattern.toStdString() << std::endl;
    }
    
    void applyDarkFrameSubtraction(EnhancedImageData& img) {
        QString darkFrame = findMatchingDarkFrame(img);
        
        if (darkFrame.isEmpty()) {
            std::cout << "  No matching dark frame found for:" << std::endl;
            std::cout << "    Exposure: " << img.exposure_time << "s" << std::endl;
            std::cout << "    Temperature: " << img.temperature << "K" << std::endl;
            std::cout << "    Binning: " << img.binning.toStdString() << std::endl;
            std::cout << "    Bayer: " << img.bayer_pattern.toStdString() << std::endl;
            std::cout << "  Using raw image" << std::endl;
            img.image = img.rawImage.clone();
            return;
        }
        
        // Load dark frame
        cv::Mat darkImage = loadDarkFrame(darkFrame);
        if (darkImage.empty()) {
            std::cout << "  Failed to load dark frame - using raw image" << std::endl;
            img.image = img.rawImage.clone();
            return;
        }
        
        // Ensure dark frame has same size
        if (darkImage.size() != img.rawImage.size()) {
            std::cout << "  Dark frame size mismatch - using raw image" << std::endl;
            img.image = img.rawImage.clone();
            return;
        }
        
        // Perform dark subtraction
        img.image = img.rawImage - darkImage;
        img.darkSubtracted = true;
        img.matchedDarkFrame = QFileInfo(darkFrame).fileName();
        
        std::cout << "  Dark subtracted: " << img.matchedDarkFrame.toStdString() << std::endl;
    }
    
    QString findMatchingDarkFrame(const EnhancedImageData& img) {
        if (m_darkDirectory.isEmpty()) {
            return QString();
        }
        
        QDir darkDir(m_darkDirectory);
        if (!darkDir.exists()) {
            return QString();
        }
        
        QStringList darkFiles = darkDir.entryList(QStringList() << "master*.fits", QDir::Files);
        
        QString bestMatch;
        
        for (const QString& darkFile : darkFiles) {
	    bestMatch = darkDir.absoluteFilePath(darkFile);
        }
        
        return bestMatch;
    }
    
    bool extractDarkMetadata(const QString& darkPath, int& exposure, int& temperature, 
                           QString& binning, QString& bayerPattern) {
        fitsfile* fptr = nullptr;
        int status = 0;
        
        QByteArray pathBytes = darkPath.toLocal8Bit();
        if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
            return false;
        }
        
        // Extract exposure (handle milliseconds)
        double exp = 0.0;
        if (fits_read_key(fptr, TDOUBLE, "EXPTIME", &exp, nullptr, &status) == 0) {
            // EXPTIME usually in seconds
            exposure = static_cast<int>(exp);
        } else if (fits_read_key(fptr, TDOUBLE, "EXPOSURE", &exp, nullptr, &status) == 0) {
            // EXPOSURE might be in milliseconds for Stellina
            if (exp > 100) {
                // Likely milliseconds, convert to seconds
                exposure = static_cast<int>(exp / 1000.0);
            } else {
                // Already in seconds
                exposure = static_cast<int>(exp);
            }
        } else {
            fits_close_file(fptr, &status);
            return false;
        }
        status = 0;
        
        // Extract temperature (handle Celsius)
        double temp = 15.0; // Default 15°C
        if (fits_read_key(fptr, TDOUBLE, "TEMP", &temp, nullptr, &status) == 0 ||
            fits_read_key(fptr, TDOUBLE, "CCD-TEMP", &temp, nullptr, &status) == 0) {
            // For Stellina, TEMP is in Celsius, convert to Kelvin
            if (temp < 100) temp += 273.15;
            temperature = static_cast<int>(temp);
        }
        status = 0;
        
        // Extract binning
        char bin[FLEN_VALUE] = "1x1";
        fits_read_key(fptr, TSTRING, "XBINNING", bin, nullptr, &status);
        binning = QString::fromLatin1(bin).trimmed().remove('\'').remove('"');
        status = 0;
        
        // Extract bayer pattern
        char bayer[FLEN_VALUE] = "RGGB";
        fits_read_key(fptr, TSTRING, "BAYERPAT", bayer, nullptr, &status);
        bayerPattern = QString::fromLatin1(bayer).trimmed().remove('\'').remove('"').toUpper();
        
        fits_close_file(fptr, &status);
        return true;
    }
    
    cv::Mat loadDarkFrame(const QString& darkPath) {
        fitsfile* fptr = nullptr;
        int status = 0;
        
        QByteArray pathBytes = darkPath.toLocal8Bit();
        if (fits_open_file(&fptr, pathBytes.data(), READONLY, &status)) {
            return cv::Mat();
        }
        
        long naxes[2];
        if (fits_get_img_size(fptr, 2, naxes, &status)) {
            fits_close_file(fptr, &status);
            return cv::Mat();
        }
        
        long totalPixels = naxes[0] * naxes[1];
        std::vector<float> pixels(totalPixels);
        
        if (fits_read_img(fptr, TFLOAT, 1, totalPixels, nullptr, pixels.data(), nullptr, &status)) {
            fits_close_file(fptr, &status);
            return cv::Mat();
        }
        
        fits_close_file(fptr, &status);
        
        return cv::Mat(naxes[1], naxes[0], CV_32F, pixels.data()).clone();
    }
    
    // Dark frame rotation methods (adapted from StellinaProcessor)
    QString createRotatedDarkFrame(const QString& originalDarkPath, const QString& fromPattern, const QString& toPattern) {
        // Generate rotated dark frame filename
        QFileInfo darkInfo(originalDarkPath);
        QString rotatedName = QString("rotated_%1_%2_to_%3.%4")
            .arg(darkInfo.baseName())
            .arg(fromPattern)
            .arg(toPattern)
            .arg(darkInfo.suffix());
        
        // Create temporary rotated dark in current directory
        QString rotatedPath = rotatedName;
        
        // Check if already exists
        if (QFile::exists(rotatedPath)) {
            return rotatedPath;
        }
        
        // Load original dark frame
        cv::Mat originalDark = loadDarkFrame(originalDarkPath);
        if (originalDark.empty()) {
            std::cout << "    Failed to load original dark for rotation" << std::endl;
            return QString();
        }
        
        // Perform rotation based on Bayer pattern transformation
        cv::Mat rotatedDark;
        if (rotateDarkForBayerPattern(originalDark, rotatedDark, fromPattern, toPattern)) {
            // Save rotated dark frame
            if (saveDarkFrame(rotatedDark, rotatedPath, originalDarkPath, fromPattern, toPattern)) {
                return rotatedPath;
            }
        }
        
        return QString();
    }
    
    bool rotateDarkForBayerPattern(const cv::Mat& inputDark, cv::Mat& outputDark, 
                                   const QString& fromPattern, const QString& toPattern) {
        if (fromPattern == toPattern) {
            outputDark = inputDark.clone();
            return true;
        }
        
        // Bayer pattern transformations with 180° rotation:
        // RGGB -> BGGR (and vice versa)
        // GRBG -> GBRG (and vice versa)
        
        if ((fromPattern == "RGGB" && toPattern == "BGGR") ||
            (fromPattern == "BGGR" && toPattern == "RGGB") ||
            (fromPattern == "GRBG" && toPattern == "GBRG") ||
            (fromPattern == "GBRG" && toPattern == "GRBG")) {
            
            std::cout << "    Rotating 180°: " << fromPattern.toStdString() 
                      << " -> " << toPattern.toStdString() << std::endl;
            return rotateDark180(inputDark, outputDark);
        }
        
        // Add other rotations as needed (90°, 270°, flips)
        std::cout << "    Unsupported Bayer rotation: " << fromPattern.toStdString() 
                  << " -> " << toPattern.toStdString() << std::endl;
        outputDark = inputDark.clone();
        return false;
    }
    
    bool rotateDark180(const cv::Mat& input, cv::Mat& output) {
        int height = input.rows;
        int width = input.cols;
        
        output = cv::Mat::zeros(height, width, input.type());
        
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int srcY = height - 1 - y;
                int srcX = width - 1 - x;
                output.at<float>(y, x) = input.at<float>(srcY, srcX);
            }
        }
        
        return true;
    }
    
    bool saveDarkFrame(const cv::Mat& darkImage, const QString& outputPath, 
                      const QString& originalPath, const QString& fromPattern, const QString& toPattern) {
        fitsfile *outputFits = nullptr;
        int status = 0;
        
        // Create output FITS file
        QByteArray pathBytes = QString("!%1").arg(outputPath).toLocal8Bit();
        if (fits_create_file(&outputFits, pathBytes.data(), &status)) {
            std::cout << "    Failed to create rotated dark frame file" << std::endl;
            return false;
        }
        
        // Create image structure
        long naxes[2] = {darkImage.cols, darkImage.rows};
        if (fits_create_img(outputFits, FLOAT_IMG, 2, naxes, &status)) {
            fits_close_file(outputFits, &status);
            return false;
        }
        
        // Convert and write image data
        std::vector<float> pixels(darkImage.rows * darkImage.cols);
        for (int y = 0; y < darkImage.rows; ++y) {
            for (int x = 0; x < darkImage.cols; ++x) {
                pixels[y * darkImage.cols + x] = darkImage.at<float>(y, x);
            }
        }
        
        if (fits_write_img(outputFits, TFLOAT, 1, pixels.size(), pixels.data(), &status)) {
            fits_close_file(outputFits, &status);
            return false;
        }
        
        // Copy header from original dark frame
        fitsfile *originalFits = nullptr;
        QByteArray originalPathBytes = originalPath.toLocal8Bit();
        if (fits_open_file(&originalFits, originalPathBytes.data(), READONLY, &status) == 0) {
            fits_copy_header(originalFits, outputFits, &status);
            fits_close_file(originalFits, &status);
            status = 0; // Reset status in case copy failed
        }
        
        // Update Bayer pattern in header
        QByteArray patternBytes = toPattern.toLocal8Bit();
        char* patternPtr = patternBytes.data();
        fits_update_key(outputFits, TSTRING, "BAYERPAT", patternPtr, "Rotated bayer pattern", &status);
        if (status != 0) {
            status = 0;
            fits_write_key(outputFits, TSTRING, "BAYERPAT", patternPtr, "Rotated bayer pattern", &status);
        }
        
        // Add processing history
        QString historyComment = QString("HISTORY Rotated dark frame from %1 to %2 pattern for Stellina")
                                    .arg(fromPattern).arg(toPattern);
        QByteArray historyBytes = historyComment.toLocal8Bit();
        fits_write_history(outputFits, historyBytes.data(), &status);
        
        QString methodComment = QString("HISTORY Rotation method: 180 degrees");
        QByteArray methodBytes = methodComment.toLocal8Bit();
        fits_write_history(outputFits, methodBytes.data(), &status);
        
        fits_close_file(outputFits, &status);
        
        if (status != 0) {
            QFile::remove(outputPath);
            return false;
        }
        
        std::cout << "    Saved rotated dark: " << QFileInfo(outputPath).fileName().toStdString() << std::endl;
        return true;
    }
    
public:
    bool stackImages() {
        if (images.empty()) {
            std::cout << "No images to stack" << std::endl;
            return false;
        }
        
        std::cout << "\n=== CONFIGURABLE STACKING ===" << std::endl;
        
        // Use first image as reference
        const EnhancedImageData& reference = images[0];
        cv::Mat stacked = cv::Mat::zeros(reference.size, CV_32F);
        cv::Mat weight_map = cv::Mat::zeros(reference.size, CV_32F);
        
        std::cout << "Reference image: " << QFileInfo(reference.filename).fileName().toStdString() << std::endl;
        std::cout << "Reference target position: (" << reference.target_pixel.x() 
                  << ", " << reference.target_pixel.y() << ")" << std::endl;
        
        // Set reference time for sidereal tracking
        if (config.mode == TrackingMode::SIDEREAL && config.reference_jd == 0.0) {
            config.reference_jd = reference.jd;
            std::cout << "Setting reference JD: " << config.reference_jd << std::endl;
        }
        
        for (size_t i = 0; i < images.size(); ++i) {
            EnhancedImageData& img = images[i];
            
            if (!img.hasValidWCS) {
                std::cout << "Skipping image " << (i+1) << " - no valid WCS" << std::endl;
                continue;
            }
            
            // Calculate transformation to align target positions
            cv::Point2f ref_center(reference.target_pixel.x(), reference.target_pixel.y());
            cv::Point2f img_center(img.target_pixel.x(), img.target_pixel.y());
            
            // Calculate field rotation correction
            double rotation_diff = 0.0;
            if (config.apply_rotation_correction) {
                rotation_diff = calculateFieldRotation(reference, img);
                img.calculated_rotation = rotation_diff;
            }
            
            // Create transformation matrix
            cv::Mat transform = calculateTransform(img_center, ref_center, rotation_diff, img.size);
            
            // Store transformation info for diagnostics
            img.translation_offset = QPointF(ref_center.x - img_center.x, ref_center.y - img_center.y);
            
            std::cout << "Image " << (i+1) << " transformation:" << std::endl;
            std::cout << "  Translation: (" << img.translation_offset.x()
                      << ", " << img.translation_offset.y() << ")" << std::endl;
            if (config.apply_rotation_correction) {
                std::cout << "  Rotation: " << rotation_diff << "°" << std::endl;
            }
            
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
        
        return saveFITSImageWithConfig(normalized);
    }
    
private:
    double calculateFieldRotation(const EnhancedImageData& ref, const EnhancedImageData& img) {
        if (!ref.hasValidWCS || !img.hasValidWCS) {
            return 0.0; // No rotation correction possible
        }
        
        // For sidereal tracking, calculate field rotation due to Earth's rotation
        if (config.mode == TrackingMode::SIDEREAL) {
            // Field rotation = hour angle difference * cos(declination)
            double time_diff_hours = (img.jd - ref.jd) * 24.0;
            double hour_angle_diff = time_diff_hours * 15.0; // degrees
            double field_rotation = hour_angle_diff * cos(config.target_dec * DEG2RAD);
            
            // Limit to reasonable values
            while (field_rotation > 180.0) field_rotation -= 360.0;
            while (field_rotation < -180.0) field_rotation += 360.0;
            
            return field_rotation;
        }
        
        // For other modes, use sky coordinate method
        QPointF cornerPixel(ref.size.width * 0.8, ref.size.height * 0.2);
        QPointF cornerSky = ref.wcsValidator->pixelToSky(cornerPixel.x(), cornerPixel.y());
        
        if (cornerSky.x() < 0) {
            return 0.0; // Invalid transformation
        }
        
        QPointF cornerInImg = img.wcsValidator->skyToPixel(cornerSky.x(), cornerSky.y());
        
        if (cornerInImg.x() < 0) {
            return 0.0; // Invalid transformation
        }
        
        // Calculate rotation based on vector differences
        QPointF refVector = cornerPixel - QPointF(ref.target_pixel);
        QPointF imgVector = cornerInImg - QPointF(img.target_pixel);
        
        double refAngle = atan2(refVector.y(), refVector.x()) * RAD2DEG;
        double imgAngle = atan2(imgVector.y(), imgVector.x()) * RAD2DEG;
        
        double rotation = refAngle - imgAngle;
        
        // Normalize to [-180, 180]
        while (rotation > 180.0) rotation -= 360.0;
        while (rotation < -180.0) rotation += 360.0;
        
        return rotation;
    }
    
    cv::Mat calculateTransform(cv::Point2f src_center, cv::Point2f dst_center, 
                              double rotation_degrees, cv::Size img_size) {
        
        // Image center for rotation
        cv::Point2f image_center(img_size.width / 2.0f, img_size.height / 2.0f);
        
        // Create rotation matrix around image center
        cv::Mat rotation_matrix = cv::getRotationMatrix2D(image_center, rotation_degrees, 1.0);
        
        // Add translation to align target centers
        cv::Point2f translation = dst_center - src_center;
        rotation_matrix.at<double>(0, 2) += translation.x;
        rotation_matrix.at<double>(1, 2) += translation.y;
        
        return rotation_matrix;
    }
    
    bool saveFITSImageWithConfig(const cv::Mat& image) {
        fitsfile *fptr = nullptr;
        int status = 0;
        
        QByteArray pathBytes = QString("!%1").arg(config.output_file).toLocal8Bit();
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
        
        // Copy WCS from reference
        if (!images.empty() && images[0].hasValidWCS) {
            const auto& reference = images[0];
            
            double centerRA, centerDec;
            if (reference.wcsValidator->getCenter(centerRA, centerDec)) {
                fits_write_key(fptr, TDOUBLE, "CRVAL1", &centerRA, "Reference RA", &status);
                fits_write_key(fptr, TDOUBLE, "CRVAL2", &centerDec, "Reference Dec", &status);
            }
            
            double crpix1 = reference.size.width / 2.0 + 1.0;
            double crpix2 = reference.size.height / 2.0 + 1.0;
            fits_write_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "Reference pixel X", &status);
            fits_write_key(fptr, TDOUBLE, "CRPIX2", &crpix2, "Reference pixel Y", &status);
            
            double pixelScale = reference.wcsValidator->getPixScale();
            double cdelt1 = -pixelScale / 3600.0;
            double cdelt2 = pixelScale / 3600.0;
            fits_write_key(fptr, TDOUBLE, "CDELT1", &cdelt1, "Pixel scale X", &status);
            fits_write_key(fptr, TDOUBLE, "CDELT2", &cdelt2, "Pixel scale Y", &status);
        }
        
        char ctype1[] = "RA---TAN";
        char ctype2[] = "DEC--TAN";
        fits_write_key(fptr, TSTRING, "CTYPE1", ctype1, "Coordinate type", &status);
        fits_write_key(fptr, TSTRING, "CTYPE2", ctype2, "Coordinate type", &status);
        
        // Add configuration info
        int nimages = images.size();
        fits_write_key(fptr, TINT, "NSTACKED", &nimages, "Number of stacked images", &status);
        
        QString methodStr;
        switch(config.mode) {
            case TrackingMode::SIDEREAL: methodStr = "SIDEREAL_TRACKING"; break;
            case TrackingMode::COMET: methodStr = "COMET_TRACKING"; break;
            case TrackingMode::FIXED: methodStr = "FIXED_REFERENCE"; break;
        }
        
        QByteArray methodBytes = methodStr.toLocal8Bit();
        char* methodPtr = methodBytes.data();
        fits_write_key(fptr, TSTRING, "STACKMET", methodPtr, "Stacking method", &status);
        
        if (config.mode == TrackingMode::SIDEREAL) {
            fits_write_key(fptr, TDOUBLE, "TARG_RA", &config.target_ra, "Target RA for tracking", &status);
            fits_write_key(fptr, TDOUBLE, "TARG_DEC", &config.target_dec, "Target Dec for tracking", &status);
            fits_write_key(fptr, TDOUBLE, "REF_JD", &config.reference_jd, "Reference Julian Date", &status);
        }
        
        fits_close_file(fptr, &status);
        
        std::cout << "Stacked image saved: " << config.output_file.toStdString() << std::endl;
        return status == 0;
    }
};

int main(int argc, char* argv[]) {
    QCoreApplication app(argc, argv);
    
    if (argc < 4) {
        std::cout << "Configurable Star/Comet Stacker with PCL WCS and Dark Subtraction" << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "  " << argv[0] << " sidereal <ra> <dec> <image_dir> [dark_dir]" << std::endl;
        std::cout << "  " << argv[0] << " comet <ephemeris.csv> <image_dir> [dark_dir]" << std::endl;
        std::cout << "  " << argv[0] << " fixed <image_dir> [dark_dir]" << std::endl;
        std::cout << std::endl;
        std::cout << "Tracking modes:" << std::endl;
        std::cout << "  sidereal: Track at sidereal rate (stars stationary, STRAIGHT trails)" << std::endl;
        std::cout << "  comet:    Track comet motion (comet stationary, star trails)" << std::endl;
        std::cout << "  fixed:    No tracking (everything moves)" << std::endl;
        std::cout << std::endl;
        std::cout << "Dark frames:" << std::endl;
        std::cout << "  Optional dark_dir for automatic dark frame subtraction" << std::endl;
        std::cout << "  Matches by exposure time, temperature, binning, and bayer pattern" << std::endl;
        std::cout << "  Supports reversed Stellina images (img-####r pattern)" << std::endl;
        std::cout << std::endl;
        std::cout << "Examples:" << std::endl;
        std::cout << "  # Sidereal tracking with dark subtraction" << std::endl;
        std::cout << "  " << argv[0] << " sidereal 136.067 79.6577 /path/to/lights /path/to/darks" << std::endl;
        std::cout << "  " << std::endl;
        std::cout << "  # Comet tracking without darks" << std::endl;
        std::cout << "  " << argv[0] << " comet ephemeris.csv /path/to/lights" << std::endl;
        return 1;
    }
    
    QString mode = QString::fromStdString(argv[1]).toLower();
    
    ConfigurableStacker stacker;
    
    // Configure based on mode
    if (mode == "sidereal" && argc >= 5) {
        double ra = std::stod(argv[2]);
        double dec = std::stod(argv[3]);
        QString image_dir = QString::fromStdString(argv[4]);
        QString dark_dir = (argc >= 6) ? QString::fromStdString(argv[5]) : QString();
        
        stacker.setConfig(TrackingMode::SIDEREAL, ra, dec, "", dark_dir);
        
        // Find FITS files
        QDir dir(image_dir);
        QStringList filters;
        filters << "*.fits" << "*.fit" << "*.FITS" << "*.FIT";
        QStringList fits_files = dir.entryList(filters, QDir::Files, QDir::Name);
        
        if (fits_files.isEmpty()) {
            std::cout << "No FITS files found in: " << image_dir.toStdString() << std::endl;
            return 1;
        }
        
        // Convert to full paths
        QStringList full_paths;
        for (const QString& file : fits_files) {
            full_paths.append(dir.absoluteFilePath(file));
        }
        
        std::cout << "Found " << fits_files.size() << " FITS files" << std::endl;
        if (!dark_dir.isEmpty()) {
            std::cout << "Dark frame directory: " << dark_dir.toStdString() << std::endl;
        }
        
        // Load and stack images
        if (stacker.loadImages(full_paths)) {
            if (stacker.stackImages()) {
                std::cout << "\n✓ Sidereal stacking completed successfully!" << std::endl;
                std::cout << "\nStars should appear as STRAIGHT trails" << std::endl;
            } else {
                std::cout << "❌ Failed to save stacked image" << std::endl;
                return 1;
            }
        } else {
            std::cout << "❌ Failed to load images with valid WCS data" << std::endl;
            return 1;
        }
        
    } else {
        std::cout << "Invalid arguments. See usage above." << std::endl;
        return 1;
    }
    
    return 0;
}

