#ifndef STAR_CATALOG_VALIDATOR_H
#define STAR_CATALOG_VALIDATOR_H

#include <QObject>
#include <QVector>
#include <QPoint>
#include <QPointF>
#include <QString>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QJsonObject>
#include <QJsonArray>
#include <memory>

#include <pcl/WCSKeywords.h>
#include <pcl/LinearTransformation.h>
#include <pcl/AstrometricMetadata.h>
// #include <pcl/StarGenerator.h>
// #include <pcl/DynamicArray.h>
#include <curl/curl.h>
#include "ImageReader.h"

// Enhanced matching structures
struct StarMatchingParameters {
    // Distance-based matching
    double maxPixelDistance = 5.0;           // Maximum pixel distance for matching
    double searchRadius = 10.0;              // Search radius in pixels
    
    // Magnitude-based matching  
    double maxMagnitudeDifference = 2.0;     // Maximum magnitude difference
    bool useMagnitudeWeighting = true;       // Weight matches by magnitude similarity
    
    // Pattern matching
    bool useTriangleMatching = true;         // Use triangle pattern matching
    int minTriangleStars = 6;               // Minimum stars for triangle matching
    double triangleTolerancePercent = 5.0;  // Triangle similarity tolerance (%)
    
    // Geometric validation
    bool useDistortionModel = true;         // Apply distortion correction
    double maxDistortionPixels = 2.0;       // Maximum allowed distortion
    
    // Quality filtering
    double minMatchConfidence = 0.7;        // Minimum confidence for valid match
    int minMatchesForValidation = 5;        // Minimum matches needed for validation
    
    // Advanced parameters
    bool useProperMotionCorrection = true;  // Apply proper motion
    double targetEpoch = 2025.5;           // Target epoch for proper motion
    bool useBayesianMatching = false;       // Use Bayesian probability matching
};

struct EnhancedStarMatch {
    int detectedIndex = -1;
    int catalogIndex = -1;
    double pixelDistance = 0.0;
    double magnitudeDifference = 0.0;
    double confidence = 0.0;                // Match confidence (0-1)
    double triangleError = 0.0;            // Triangle pattern error
    QVector<int> supportingMatches;        // Indices of supporting matches
    bool isGeometricallyValid = false;     // Passed geometric validation
    bool isPhotometricallyValid = false;   // Passed photometric validation
    
    // Quality metrics
    double snr = 0.0;                      // Signal-to-noise ratio
    double fwhm = 0.0;                     // Full width half maximum
    double ellipticity = 0.0;              // Star ellipticity
};

struct TrianglePattern {
    QVector<int> starIndices;              // Three star indices
    double side1, side2, side3;            // Triangle side lengths
    double angle1, angle2, angle3;         // Triangle angles
    double area;                           // Triangle area
    QPointF centroid;                      // Triangle centroid
    
    // Normalized invariants (scale/rotation independent)
    double ratio12, ratio13, ratio23;      // Side ratios
    double normalizedArea;                 // Area normalized by perimeter^2
};

struct StarMatch {
    int detectedIndex = -1;     // Index in detected stars array
    int catalogIndex = -1;      // Index in catalog stars array
    double distance = 0.0;      // Distance in pixels
    double magnitudeDiff = 0.0; // Magnitude difference
    bool isGoodMatch = false;   // Whether this is considered a good match
    
    StarMatch() = default;
    StarMatch(int detIdx, int catIdx, double dist, double magDiff = 0.0)
        : detectedIndex(detIdx), catalogIndex(catIdx), distance(dist), 
          magnitudeDiff(magDiff), isGoodMatch(true) {}
};

struct CatalogStar {
    QString id;             // Star identifier (e.g., HD123456, HIP67890)
    double ra = 0.0;        // Right ascension (degrees)
    double dec = 0.0;       // Declination (degrees)
    double magnitude = 0.0; // Visual magnitude
    QString spectralType;   // Spectral type (if available)
    QPointF pixelPos;       // Calculated pixel position
    bool isValid = false;

    // NEW: Add optional enhanced photometry fields
    double magBP = 0.0;
    double magRP = 0.0;
    double colorBP_RP = 0.0;
    bool hasExtendedData = false;

    CatalogStar() = default;
    CatalogStar(const QString& starId, double rightAsc, double declin, double mag)
        : id(starId), ra(rightAsc), dec(declin), magnitude(mag), isValid(true) {}
};

struct ValidationResult {
    QVector<CatalogStar> catalogStars;
    QVector<StarMatch> matches;
    QVector<int> unmatchedDetected;   // Indices of unmatched detected stars
    QVector<int> unmatchedCatalog;    // Indices of unmatched catalog stars
    
    // Statistics
    int totalDetected = 0;
    int totalCatalog = 0;
    int totalMatches = 0;
    double matchPercentage = 0.0;
    double averagePositionError = 0.0;
    double rmsPositionError = 0.0;
    QString summary;
    bool isValid = false;
    
    void clear() {
        catalogStars.clear();
        matches.clear();
        unmatchedDetected.clear();
        unmatchedCatalog.clear();
        totalDetected = totalCatalog = totalMatches = 0;
        matchPercentage = averagePositionError = rmsPositionError = 0.0;
        summary.clear();
        isValid = false;
    }
};

// Enhanced validation result
struct EnhancedValidationResult : public ValidationResult {
    QVector<EnhancedStarMatch> enhancedMatches;
    QVector<TrianglePattern> detectedTriangles;
    QVector<TrianglePattern> catalogTriangles;
    QVector<QPair<int, int>> triangleMatches;  // Matched triangle pairs
    
    // Advanced statistics
    double geometricRMS = 0.0;             // RMS of geometric residuals
    double photometricRMS = 0.0;           // RMS of photometric residuals
    double astrometricAccuracy = 0.0;      // Estimated astrometric accuracy
    double matchingConfidence = 0.0;       // Overall matching confidence
    double scaleError = 0.0;               // Scale error percentage
    double rotationError = 0.0;            // Rotation error in degrees
    
    // Distortion analysis
    QVector<QPointF> residualVectors;      // Residual vectors for each match
    QVector<double> radialDistortions;     // Radial distortion measurements
    
    QString enhancedSummary;               // Detailed analysis summary
};

class EnhancedStarMatcher {
public:
    explicit EnhancedStarMatcher(const StarMatchingParameters& params = StarMatchingParameters());
    
    // Main matching method using PixInsight algorithms
    EnhancedValidationResult matchStarsAdvanced(
        const QVector<QPoint>& detectedStars,
        const QVector<float>& detectedMagnitudes,
        const QVector<CatalogStar>& catalogStars,
        const pcl::AstrometricMetadata& astrometry);
    
    // Individual matching algorithms
    QVector<EnhancedStarMatch> performInitialMatching(
        const QVector<QPoint>& detected,
        const QVector<CatalogStar>& catalog,
        const pcl::AstrometricMetadata& astrometry);
        
    QVector<EnhancedStarMatch> performTriangleMatching(
        const QVector<QPoint>& detected,
        const QVector<CatalogStar>& catalog);
        
    QVector<EnhancedStarMatch> performGeometricValidation(
        const QVector<EnhancedStarMatch>& initialMatches,
        const QVector<QPoint>& detected,
        const QVector<CatalogStar>& catalog);
    
    // Pattern analysis
    QVector<TrianglePattern> generateTrianglePatterns(
        const QVector<QPoint>& stars, bool isPixelCoords = true);
    QVector<TrianglePattern> generateTrianglePatterns(
        const QVector<CatalogStar>& stars, bool isPixelCoords = false);
        
    QPair<QVector<QPair<int, int>>, double> matchTrianglePatterns(
        const QVector<TrianglePattern>& pattern1,
        const QVector<TrianglePattern>& pattern2);
    
    // Distortion analysis
    void analyzeDistortions(EnhancedValidationResult& result,
                           const QVector<QPoint>& detected,
                           const QVector<CatalogStar>& catalog);
    
    // Quality assessment
    double calculateMatchConfidence(const EnhancedStarMatch& match,
                                   const QVector<EnhancedStarMatch>& allMatches);
    
    void filterLowQualityMatches(QVector<EnhancedStarMatch>& matches);
    
    // Utility methods
    void setParameters(const StarMatchingParameters& params) { m_params = params; }
    StarMatchingParameters getParameters() const { return m_params; }
    
private:
    StarMatchingParameters m_params;
    
    // Internal methods
    double calculateTriangleSimilarity(const TrianglePattern& t1, const TrianglePattern& t2);
    bool isTriangleValid(const TrianglePattern& triangle);
    QPointF calculateResidualVector(const EnhancedStarMatch& match,
                                   const QVector<QPoint>& detected,
                                   const QVector<CatalogStar>& catalog);
    double estimateMagnitudeFromFlux(double flux, const QVector<float>& allFluxes);
    
    // Statistical analysis
    void calculateAdvancedStatistics(EnhancedValidationResult& result);
    double calculateGeometricRMS(const QVector<EnhancedStarMatch>& matches,
                                const QVector<QPoint>& detected,
                                const QVector<CatalogStar>& catalog);
};
    
struct WCSData {
    double crval1 = 0.0;    // Reference RA (degrees)
    double crval2 = 0.0;    // Reference Dec (degrees)
    double crpix1 = 0.0;    // Reference pixel X
    double crpix2 = 0.0;    // Reference pixel Y
    double cd11 = 0.0;      // CD matrix elements
    double cd12 = 0.0;
    double cd21 = 0.0;
    double cd22 = 0.0;
    double pixscale = 0.0;  // Pixel scale (arcsec/pixel)
    double orientation = 0.0; // Position angle (degrees)
    int width = 0;          // Image width
    int height = 0;         // Image height
    bool isValid = false;
    
    // Alternative: CDELT/PC matrix format (less common but supported)
    double cdelt1 = 0.0;
    double cdelt2 = 0.0;
    double pc11 = 1.0;
    double pc12 = 0.0;
    double pc21 = 0.0;
    double pc22 = 1.0;
    
    void clear() {
        crval1 = crval2 = crpix1 = crpix2 = 0.0;
        cd11 = cd12 = cd21 = cd22 = 0.0;
        pixscale = orientation = 0.0;
        width = height = 0;
        isValid = false;
        cdelt1 = cdelt2 = 0.0;
        pc11 = pc22 = 1.0;
        pc12 = pc21 = 0.0;
    }
      
    void reset() {
        isValid = false;
        crval1 = crval2 = crpix1 = crpix2 = 0.0;
        cd11 = cd12 = cd21 = cd22 = 0.0;
        pixscale = orientation = 0.0;
        width = height = 0;
    }
    WCSData() = default;

};

class StarCatalogValidator : public QObject
{
    Q_OBJECT
    
public:
    enum CatalogSource {
        Hipparcos,      // Bright stars (mag < 12)
        Tycho2,         // Fainter stars (mag < 12)
        Gaia,           // Very comprehensive (millions of stars)
        Custom          // User-provided catalog
    };
    
    enum ValidationMode {
        Strict,         // Require very close matches
        Normal,         // Standard matching tolerance
        Loose          // More permissive matching
    };
    
    explicit StarCatalogValidator(QObject* parent = nullptr);
    ~StarCatalogValidator();
    
    // Enhanced validation method
    EnhancedValidationResult validateStarsAdvanced(
        const QVector<QPoint>& detectedStars,
        const QVector<float>& starMagnitudes = QVector<float>(),
        const StarMatchingParameters& params = StarMatchingParameters());
    
    // Set advanced matching parameters
    void setMatchingParameters(const StarMatchingParameters& params);
    StarMatchingParameters getMatchingParameters() const;
    
    // Pattern-based catalog search
    void queryCatalogWithPatterns(double centerRA, double centerDec, 
                                 double radiusDegrees,
                                 const QVector<QPoint>& referencePattern);
    
    // Distortion modeling
    bool calibrateDistortionModel(const EnhancedValidationResult& validationResult);
    QPointF applyDistortionCorrection(const QPointF& pixelPos);
    // Configuration
    void setCatalogSource(CatalogSource source);
    void setValidationMode(ValidationMode mode);
    void setMatchingTolerance(double pixelTolerance, double magnitudeTolerance = 2.0);
  //    void setMagnitudeLimit(double faintestMagnitude);
    
    // WCS handling
    bool setWCSData(const WCSData& wcs);
    void setWCSFromMetadata(const QStringList& metadata);
    WCSData getWCSData() const;
    bool hasValidWCS() const;
    
    // Coordinate transformations
    bool setWCSFromPCLKeywords(const pcl::FITSKeywordArray& keywords);
    bool setWCSFromImageMetadata(const ImageData& imageData);  // NEW: Direct from ImageData
    bool setWCSFromSolver(pcl::FITSKeywordArray keywords, int width, int height);
    
    // PCL-based coordinate transformations
    QPointF skyToPixel(double ra, double dec) const;
    QPointF pixelToSky(double x, double y) const;
    
    // Catalog operations
    void queryCatalog(double centerRA, double centerDec, double radiusDegrees);
    void loadCustomCatalog(const QString& filePath);
    void loadCustomCatalog(const QVector<CatalogStar>& stars);
    void addBrightStarsFromDatabase(double centerRA, double centerDec, double radiusDegrees);
    void queryGaiaDR3(double centerRA, double centerDec, double radiusDegrees);
    void queryGaiaWithSpectra(double centerRA, double centerDec, double radiusDegrees);  
    void findBrightGaiaStars(double centerRA, double centerDec, double radiusDegrees, int count = 20);
    void initializeGaiaDR3();
  
    // Validation
    ValidationResult validateStars(const QVector<QPoint>& detectedStars, 
                                   const QVector<float>& starRadii = QVector<float>());
    
    // Results access
    ValidationResult getLastValidation() const { return m_lastValidation; }
    QVector<CatalogStar> getCatalogStars();
    void showCatalogStats() const;
    
    // Utility
    QString getValidationSummary() const;
    QJsonObject exportValidationResults() const;
    void clearResults();
    void showWCSInfoFromPCL();
    void extractAndDisplayCDMatrix();
    void calculateAndDisplayFOV();
    void performWCSQualityCheck()    ;
    void displayPCLSpecificInfo();
    void showBasicWCSInfoFromPCL();
  double getPixScale()
  {
    return m_astrometricMetadata.Resolution() * 3600.0; // arcsec/px
  }
    bool getCenter(double &x, double &y)
    {
      pcl::DPoint centerCoords;
      if (m_astrometricMetadata.ImageCenterToCelestial(centerCoords)) {
	x = centerCoords.x;
	y = centerCoords.y;
	return true;
      }
      return false;
    }
    int getWidth() { return m_astrometricMetadata.Width(); }
    int getHeight() { return m_astrometricMetadata.Height(); }
    void setMetadata(pcl::AstrometricMetadata rslt) { m_astrometricMetadata = rslt; }
    void testAstrometricMetadata() const;
  
signals:
    void catalogQueryStarted();
    void catalogQueryFinished(bool success, const QString& message);
    void validationCompleted(const ValidationResult& result);
    void errorSignal(const QString& message);
    
private slots:
  //    void onCatalogQueryFinished();
    void onNetworkError(QNetworkReply::NetworkError error);
    
private:

    std::unique_ptr<EnhancedStarMatcher> m_enhancedMatcher;
    StarMatchingParameters m_matchingParams;
    
    // Distortion model parameters
    bool m_hasDistortionModel = false;
    double m_radialDistortionK1 = 0.0;
    double m_radialDistortionK2 = 0.0;
    double m_tangentialDistortionP1 = 0.0;
    double m_tangentialDistortionP2 = 0.0;
    QPointF m_distortionCenter;

    // Replace complex WCS handling with PCL's high-level class
    pcl::AstrometricMetadata m_astrometricMetadata;
  //    bool m_hasAstrometricData = false;
    
    // Helper methods
    // Configuration
    ValidationMode m_validationMode;
    double m_pixelTolerance;
    double m_magnitudeTolerance;
  //    double m_magnitudeLimit;
    
    // Data
  //    WCSData m_wcsData;
    QVector<CatalogStar> m_catalogStars;
    ValidationResult m_lastValidation;
    
    // Network
    CURL* m_curl;
    QByteArray m_catData;
    
    // Internal methods
    void initializeTolerances();
    QString buildCatalogQuery(double centerRA, double centerDec, double radiusDegrees);
    void parseCatalogResponse(const QByteArray& data);
    void parseHipparcosData(const QJsonArray& stars);
    void parseGaiaData(const QJsonArray& stars);
    double parseCoordinate(const QString& coordStr, bool isRA) const;
    
    ValidationResult performMatching(const QVector<QPoint>& detectedStars, 
                                     const QVector<float>& starRadii);
    double calculateDistance(const QPoint& detected, const QPointF& catalog) const;
    bool isGoodMatch(double distance, double magnitudeDiff) const;
    void calculateStatistics(ValidationResult& result) const;
    void queryGaiaCatalog(double centerRA, double centerDec, double radiusDegrees);
  
    // WCS transformation helpers
    void updateWCSMatrix();
    void testPCLWCS() const;
    bool m_wcsMatrixValid;
    double m_transformMatrix[4]; // 2x2 transformation matrix
    double m_det; // Determinant for inverse transformation
    // Use PCL's native WCS objects
};

#endif // STAR_CATALOG_VALIDATOR_H
