#ifndef GAIA_GDR3_CATALOG_H
#define GAIA_GDR3_CATALOG_H

#include <QDebug>
#include <QString>
#include <QVector>
#include <QTime>
#include <QFileInfo>
#include <QMutex>
#include <memory>

// Forward declare PCL classes to avoid including in header
namespace pcl {
    class GaiaDatabaseFile;
    struct GaiaSearchData;
}

class GaiaGDR3Catalog
{
public:
    struct Star {
        double ra, dec;           // Coordinates (degrees)
        double magnitude;         // G magnitude
        double magBP, magRP;      // BP and RP magnitudes
        double parallax;          // Parallax (mas)
        double pmRA, pmDec;       // Proper motion (mas/year)
        QString sourceId;         // Gaia source ID
        QString spectralClass;    // Derived from BP-RP color
        uint32_t flags;          // Quality flags
        bool hasSpectrum;        // Has BP/RP spectrum data
        bool isValid = true;
        
        Star() : ra(0), dec(0), magnitude(0), magBP(0), magRP(0),
                 parallax(0), pmRA(0), pmDec(0), flags(0), hasSpectrum(false) {}
                 
        Star(const QString& id, double rightAsc, double declin, double mag)
            : sourceId(id), ra(rightAsc), dec(declin), magnitude(mag),
              magBP(0), magRP(0), parallax(0), pmRA(0), pmDec(0), 
              flags(0), hasSpectrum(false), isValid(true) {}
    };
    
    struct SearchParameters {
        double centerRA = 0.0;
        double centerDec = 0.0;
        double radiusDegrees = 1.0;
        double maxMagnitude = 20.0;
        double minMagnitude = -2.0;
        int maxResults = 10000;
        bool requireSpectrum = false;
        bool useProperMotion = false;    // Apply proper motion to current epoch
        double epochYear = 2025.5;       // Target epoch for proper motion
        uint32_t requiredFlags = 0;      // Quality flag requirements
        uint32_t exclusionFlags = 0;     // Quality flags to exclude
        
        SearchParameters() = default;
        SearchParameters(double ra, double dec, double radius, double magLimit = 20.0)
            : centerRA(ra), centerDec(dec), radiusDegrees(radius), maxMagnitude(magLimit) {}
    };
    
private:
    static QString s_catalogPath;
    static std::unique_ptr<pcl::GaiaDatabaseFile> s_database;
    static QMutex s_databaseMutex;
    static bool s_isInitialized;
    
    // Internal methods
    static bool initializeDatabase();
    static void cleanupDatabase();
    static QString deriveSpectralClass(double bpMag, double rpMag, double gMag);
    static void applyProperMotion(Star& star, double targetEpoch);
    
public:
    // Static interface methods
    static void setCatalogPath(const QString& path);
    static bool isAvailable();
    static QString getCatalogInfo();
    
    // Main query method
    static QVector<Star> queryRegion(const SearchParameters& params);
    
    // Convenience methods
    static QVector<Star> queryRegion(double centerRA, double centerDec, 
                                   double radiusDegrees, double maxMagnitude = 20.0);
    
    static QVector<Star> findBrightestStars(double centerRA, double centerDec,
                                          double radiusDegrees, int count = 50);
    
    static QVector<Star> findStarsWithSpectra(double centerRA, double centerDec,
                                            double radiusDegrees, double maxMagnitude = 15.0);
    
    // Utility methods
    static void testPerformance(double centerRA, double centerDec, double radius);
    static void printCatalogStatistics();
    static bool validateDatabase();
    
    // Quality flag helpers
    static bool hasGoodAstrometry(const Star& star);
    static bool hasGoodPhotometry(const Star& star);
    static bool isHighQuality(const Star& star);
    
    // Coordinate utilities
    static double calculateAngularSeparation(double ra1, double dec1, double ra2, double dec2);
    static QPair<double, double> galacticCoordinates(double ra, double dec);
    
    // Cleanup (call at application exit)
    static void shutdown();
};

#endif // GAIA_GDR3_CATALOG_H