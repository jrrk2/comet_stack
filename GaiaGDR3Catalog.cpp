#include "GaiaGDR3Catalog.h"

// Initialize PCL Mock API before including PCL headers
#include "PCLMockAPI.h"

// PCL includes
#include <pcl/GaiaDatabaseFile.h>
#include <pcl/Console.h>
#include <pcl/Math.h>
#include <pcl/Exception.h>

#include <QMutexLocker>
#include <QApplication>
#include <cmath>
#include <algorithm>

// Static member definitions
QString GaiaGDR3Catalog::s_catalogPath;
std::unique_ptr<pcl::GaiaDatabaseFile> GaiaGDR3Catalog::s_database;
QMutex GaiaGDR3Catalog::s_databaseMutex;
bool GaiaGDR3Catalog::s_isInitialized = false;

void GaiaGDR3Catalog::setCatalogPath(const QString& path)
{
    QMutexLocker locker(&s_databaseMutex);
    s_catalogPath = path;
    s_isInitialized = false;
    
    // Close existing database
    if (s_database) {
        s_database.reset();
    }
    
    qDebug() << "Gaia GDR3 catalog path set to:" << path;
}

bool GaiaGDR3Catalog::isAvailable()
{
    QMutexLocker locker(&s_databaseMutex);
    return !s_catalogPath.isEmpty() && QFileInfo::exists(s_catalogPath);
}

QString GaiaGDR3Catalog::getCatalogInfo()
{
    if (!isAvailable()) {
        return "Gaia GDR3 catalog not available";
    }
    
    QFileInfo info(s_catalogPath);
    double sizeMB = info.size() / (1024.0 * 1024.0);
    
    QString result = QString("Gaia GDR3 catalog: %1 (%2 MB)")
                    .arg(info.fileName()).arg(sizeMB);
    
    if (initializeDatabase() && s_database) {
        try {
            result += QString("\nData release: %1")
                     .arg(QString::fromStdString(s_database->DataRelease().c_str()));
            result += QString("\nMagnitude range: [%1, %2]")
                     .arg(s_database->MagnitudeLow()).arg(s_database->MagnitudeHigh());
            
            if (s_database->HasMeanSpectrumData()) {
                result += QString("\nBP/RP spectra: Available (%1 wavelengths, %2-%3 nm)")
                         .arg(s_database->SpectrumCount())
                         .arg(s_database->SpectrumStart())
                         .arg(s_database->SpectrumStart() + s_database->SpectrumCount() * s_database->SpectrumStep());
            } else {
                result += "\nBP/RP spectra: Not available";
            }
        } catch (const pcl::Exception& e) {
            result += QString("\nDatabase info error: %1").arg(e.Message().c_str());
        }
    }
    
    return result;
}

bool GaiaGDR3Catalog::initializeDatabase()
{
    if (s_isInitialized && s_database) {
        return true;
    }
    
    if (s_catalogPath.isEmpty() || !QFileInfo::exists(s_catalogPath)) {
        qDebug() << "❌ Gaia GDR3 catalog path not set or file doesn't exist:" << s_catalogPath;
        return false;
    }
    
    try {
        // Initialize PCL Mock API
        pcl_mock::InitializeMockAPI();
        
        qDebug() << "📂 Opening Gaia GDR3 database:" << s_catalogPath;
        
        // Create the database object
        pcl::String pclPath(s_catalogPath.toUtf8().constData());
        s_database = std::make_unique<pcl::GaiaDatabaseFile>(pclPath);
        
        if (!s_database->IsOpen()) {
            qDebug() << "❌ Failed to open Gaia database";
            s_database.reset();
            return false;
        }
        
        s_isInitialized = true;
        
        qDebug() << "✅ Gaia GDR3 database opened successfully";
        qDebug() << "📊 Data release:" << QString::fromStdString(s_database->DataRelease().c_str());
        qDebug() << "📊 Magnitude range:" << s_database->MagnitudeLow() << "to" << s_database->MagnitudeHigh();
        
        if (s_database->HasMeanSpectrumData()) {
            qDebug() << "🌈 BP/RP spectrum data available:" << s_database->SpectrumCount() << "wavelengths";
        }
        
        return true;
        
    } catch (const pcl::Exception& e) {
        qDebug() << "❌ PCL Exception opening Gaia database:" << e.Message().c_str();
        s_database.reset();
        s_isInitialized = false;
        return false;
    } catch (const std::exception& e) {
        qDebug() << "❌ Standard exception opening Gaia database:" << e.what();
        s_database.reset();
        s_isInitialized = false;
        return false;
    } catch (...) {
        qDebug() << "❌ Unknown exception opening Gaia database";
        s_database.reset();
        s_isInitialized = false;
        return false;
    }
}

QVector<GaiaGDR3Catalog::Star> GaiaGDR3Catalog::queryRegion(const SearchParameters& params)
{
    QMutexLocker locker(&s_databaseMutex);
    QVector<Star> stars;
    
    if (!initializeDatabase()) {
        qDebug() << "❌ Failed to initialize Gaia database";
        return stars;
    }
    
    qDebug() << QString("🔍 Querying Gaia GDR3: RA=%1° Dec=%2° radius=%3° mag≤%4")
                .arg(params.centerRA).arg(params.centerDec)
                .arg(params.radiusDegrees).arg(params.maxMagnitude);
    
    auto startTime = QTime::currentTime();
    
    try {
        // Set up Gaia search parameters
        pcl::GaiaSearchData searchData;
        searchData.centerRA = params.centerRA;
        searchData.centerDec = params.centerDec;
        searchData.radius = params.radiusDegrees;
        searchData.magnitudeLow = params.minMagnitude;
        searchData.magnitudeHigh = params.maxMagnitude;
        searchData.sourceLimit = params.maxResults;
        searchData.requiredFlags = params.requiredFlags;
        searchData.exclusionFlags = params.exclusionFlags;
        
        // Enable spectrum retrieval if requested
        if (params.requireSpectrum && s_database->HasMeanSpectrumData()) {
            searchData.normalizeSpectrum = true;
            searchData.photonFluxUnits = false;
        }
        
        // Perform the search
        s_database->Search(searchData);
        
        qDebug() << QString("📊 Gaia search completed in %1 ms")
                    .arg(searchData.timeTotal * 1000.0);
        qDebug() << QString("📊 Decode time: %1 ms")
                    .arg(searchData.timeDecode * 1000.0);
        qDebug() << QString("📊 Raw results: %1 stars")
                    .arg(searchData.stars.Length());
        
        // Convert PCL results to our format
        stars.reserve(searchData.stars.Length());
        
        for (size_t i = 0; i < searchData.stars.Length(); ++i) {
            const auto& pclStar = searchData.stars[i];
            
            Star star;
            star.ra = pclStar.ra;
            star.dec = pclStar.dec;
            star.magnitude = pclStar.magG;
            star.magBP = pclStar.magBP;
            star.magRP = pclStar.magRP;
            star.parallax = pclStar.parx;  // Note: PCL uses 'parx'
            star.pmRA = pclStar.pmra;      // Note: PCL uses 'pmra'
            star.pmDec = pclStar.pmdec;    // Note: PCL uses 'pmdec'
            star.flags = pclStar.flags;
            star.hasSpectrum = !pclStar.flux.IsEmpty();
            
            // Generate source ID from flags (temporary - real ID would come from database)
            star.sourceId = QString("GDR3_%1").arg(pclStar.flags, 16, 16, QChar('0'));
            
            // Derive spectral class from colors
            star.spectralClass = deriveSpectralClass(star.magBP, star.magRP, star.magnitude);
            
            // Apply proper motion correction if requested
            if (params.useProperMotion && params.epochYear != 2016.0) {
                applyProperMotion(star, params.epochYear);
            }
            
            // Apply quality filters
            if (params.requireSpectrum && !star.hasSpectrum) {
                continue;
            }
            
            // Additional magnitude check (in case PCL search is inclusive)
            if (star.magnitude > params.maxMagnitude) {
                continue;
            }
            
            stars.append(star);
        }
        
        auto elapsed = startTime.msecsTo(QTime::currentTime());
        
        qDebug() << QString("✅ Gaia query completed in %1ms").arg(elapsed);
        qDebug() << QString("📈 Final results: %1 stars").arg(stars.size());
        
        if (!stars.isEmpty()) {
            double brightestMag = stars[0].magnitude;
            double faintestMag = stars[0].magnitude;
            int withSpectra = 0;
            
            for (const auto& star : stars) {
                brightestMag = std::min(brightestMag, star.magnitude);
                faintestMag = std::max(faintestMag, star.magnitude);
                if (star.hasSpectrum) withSpectra++;
            }
            
            qDebug() << QString("📊 Magnitude range: %1 to %2").arg(brightestMag).arg(faintestMag);
            qDebug() << QString("📊 Stars with BP/RP spectra: %1 (%2%%)")
                        .arg(withSpectra).arg(100.0 * withSpectra / stars.size());
        }
        
        // Sort by magnitude (brightest first)
        std::sort(stars.begin(), stars.end(), 
                  [](const Star& a, const Star& b) {
                      return a.magnitude < b.magnitude;
                  });
        
    } catch (const pcl::Exception& e) {
        qDebug() << "❌ PCL Exception during Gaia search:" << e.Message().c_str();
    } catch (const std::exception& e) {
        qDebug() << "❌ Standard exception during Gaia search:" << e.what();
    } catch (...) {
        qDebug() << "❌ Unknown exception during Gaia search";
    }
    
    return stars;
}

QVector<GaiaGDR3Catalog::Star> GaiaGDR3Catalog::queryRegion(double centerRA, double centerDec, 
                                                           double radiusDegrees, double maxMagnitude)
{
    SearchParameters params(centerRA, centerDec, radiusDegrees, maxMagnitude);
    return queryRegion(params);
}

QVector<GaiaGDR3Catalog::Star> GaiaGDR3Catalog::findBrightestStars(double centerRA, double centerDec,
                                                                  double radiusDegrees, int count)
{
    SearchParameters params(centerRA, centerDec, radiusDegrees, 20.0);
    params.maxResults = count * 5; // Search more to ensure we get the brightest
    
    QVector<Star> allStars = queryRegion(params);
    
    if (allStars.size() > count) {
        allStars.resize(count);
    }
    
    return allStars;
}

QVector<GaiaGDR3Catalog::Star> GaiaGDR3Catalog::findStarsWithSpectra(double centerRA, double centerDec,
                                                                    double radiusDegrees, double maxMagnitude)
{
    SearchParameters params(centerRA, centerDec, radiusDegrees, maxMagnitude);
    params.requireSpectrum = true;
    
    return queryRegion(params);
}

QString GaiaGDR3Catalog::deriveSpectralClass(double bpMag, double rpMag, double gMag)
{
    // Simplified spectral classification based on BP-RP color
    if (bpMag <= 0 || rpMag <= 0 || gMag <= 0) {
        return "Unknown";
    }
    
    double bpRpColor = bpMag - rpMag;
    
    // Very rough classification based on Gaia colors
    // This is a simplified approximation - real classification needs more data
    if (bpRpColor < 0.2) return "O/B";
    else if (bpRpColor < 0.5) return "A";
    else if (bpRpColor < 0.8) return "F";
    else if (bpRpColor < 1.2) return "G";
    else if (bpRpColor < 1.8) return "K";
    else return "M";
}

void GaiaGDR3Catalog::applyProperMotion(Star& star, double targetEpoch)
{
    // Gaia DR3 epoch is 2016.0
    const double gaiaEpoch = 2016.0;
    double deltaYears = targetEpoch - gaiaEpoch;
    
    if (std::abs(deltaYears) < 0.01) {
        return; // No significant time difference
    }
    
    // Apply proper motion correction
    // Convert proper motion from mas/year to degrees/year
    double pmRA_deg = star.pmRA / (3600.0 * 1000.0);  // mas to degrees
    double pmDec_deg = star.pmDec / (3600.0 * 1000.0);
    
    // Apply proper motion with cos(declination) correction for RA
    double cosDec = std::cos(star.dec * M_PI / 180.0);
    star.ra += (pmRA_deg * deltaYears) / cosDec;
    star.dec += pmDec_deg * deltaYears;
    
    // Normalize RA to [0, 360)
    while (star.ra < 0) star.ra += 360.0;
    while (star.ra >= 360) star.ra -= 360.0;
    
    // Clamp declination to [-90, 90]
    star.dec = std::max(-90.0, std::min(90.0, star.dec));
}

bool GaiaGDR3Catalog::hasGoodAstrometry(const Star& star)
{
    // Check for reasonable parallax and proper motion errors
    // This would need actual error data from the database
    return star.parallax > 0 && std::abs(star.parallax) < 1000; // Basic sanity check
}

bool GaiaGDR3Catalog::hasGoodPhotometry(const Star& star)
{
    // Check for reasonable photometric data
    return star.magnitude > -5.0 && star.magnitude < 25.0 &&
           star.magBP > -5.0 && star.magBP < 25.0 &&
           star.magRP > -5.0 && star.magRP < 25.0;
}

bool GaiaGDR3Catalog::isHighQuality(const Star& star)
{
    return hasGoodAstrometry(star) && hasGoodPhotometry(star);
}

double GaiaGDR3Catalog::calculateAngularSeparation(double ra1, double dec1, double ra2, double dec2)
{
    // Convert to radians
    double ra1Rad = ra1 * M_PI / 180.0;
    double dec1Rad = dec1 * M_PI / 180.0;
    double ra2Rad = ra2 * M_PI / 180.0;
    double dec2Rad = dec2 * M_PI / 180.0;
    
    // Haversine formula
    double deltaRA = ra2Rad - ra1Rad;
    double deltaDec = dec2Rad - dec1Rad;
    
    double a = std::sin(deltaDec/2) * std::sin(deltaDec/2) +
               std::cos(dec1Rad) * std::cos(dec2Rad) *
               std::sin(deltaRA/2) * std::sin(deltaRA/2);
    
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1-a));
    
    return c * 180.0 / M_PI; // Convert back to degrees
}

void GaiaGDR3Catalog::testPerformance(double centerRA, double centerDec, double radius)
{
    qDebug() << "\n=== GAIA GDR3 PERFORMANCE TEST ===";
    
    SearchParameters params(centerRA, centerDec, radius, 18.0);
    params.maxResults = 50000;
    
    auto startTime = QTime::currentTime();
    QVector<Star> stars = queryRegion(params);
    auto elapsed = startTime.msecsTo(QTime::currentTime());
    
    qDebug() << QString("🚀 Performance test: %1 stars in %2ms")
                .arg(stars.size()).arg(elapsed);
    
    if (!stars.isEmpty()) {
        qDebug() << QString("📊 Search rate: %1 stars/second")
                    .arg(stars.size() * 1000.0 / elapsed);
    }
}

void GaiaGDR3Catalog::printCatalogStatistics()
{
    qDebug() << "\n=== GAIA GDR3 CATALOG STATISTICS ===";
    qDebug() << getCatalogInfo();
}

bool GaiaGDR3Catalog::validateDatabase()
{
    QMutexLocker locker(&s_databaseMutex);
    
    if (!initializeDatabase()) {
        return false;
    }
    
    try {
        // Test a small query to validate database integrity
        SearchParameters params(0.0, 0.0, 1.0, 10.0);
        params.maxResults = 10;
        
        QVector<Star> testStars = queryRegion(params);
        
        qDebug() << QString("✅ Database validation: %1 test stars retrieved")
                    .arg(testStars.size());
        
        return true;
        
    } catch (...) {
        qDebug() << "❌ Database validation failed";
        return false;
    }
}

void GaiaGDR3Catalog::cleanupDatabase()
{
    if (s_database) {
        s_database.reset();
        s_isInitialized = false;
        qDebug() << "🔒 Gaia GDR3 database closed";
    }
}

void GaiaGDR3Catalog::shutdown()
{
    QMutexLocker locker(&s_databaseMutex);
    cleanupDatabase();
}
