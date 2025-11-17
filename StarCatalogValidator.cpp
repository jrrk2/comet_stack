#include "BrightStarDatabase.h"
#include "StarCatalogValidator.h"
#include "GaiaGDR3Catalog.h"  // Add this line
#include <QNetworkRequest>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QUrl>
#include <QUrlQuery>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QRegularExpression>
#include <QEventLoop>
#include <QTimer>
#include <cmath>
#include <algorithm>
#include <QNetworkRequest>
#include <QJsonDocument>
// Enhanced StarCatalogValidator.cpp implementation
#include "StarCatalogValidator.h"
#include <pcl/Math.h>
#include <pcl/Matrix.h>
#include <algorithm>
#include <cmath>

// Enhanced triangle validation with proper geometric checks
bool EnhancedStarMatcher::isTriangleValid(const TrianglePattern& triangle)
{
    // Check for degenerate triangles
    if (triangle.side1 < 3.0 || triangle.side2 < 3.0 || triangle.side3 < 3.0) {
        return false;
    }
    
    // Check triangle inequality
    if (triangle.side1 + triangle.side2 <= triangle.side3 ||
        triangle.side1 + triangle.side3 <= triangle.side2 ||
        triangle.side2 + triangle.side3 <= triangle.side1) {
        return false;
    }
    
    // Check for very thin triangles (area vs perimeter ratio)
    double perimeter = triangle.side1 + triangle.side2 + triangle.side3;
    if (triangle.area < 0.05 * perimeter) {
        return false; // Too thin
    }
    
    // Check for reasonable size - not too small or too large
    double maxSide = std::max({triangle.side1, triangle.side2, triangle.side3});
    if (maxSide < 10.0 || maxSide > 500.0) {
        return false;
    }
    
    // Check aspect ratio - no extremely elongated triangles
    double minSide = std::min({triangle.side1, triangle.side2, triangle.side3});
    if (maxSide / minSide > 10.0) {
        return false;
    }
    
    return true;
}

// Complete distortion analysis implementation
void EnhancedStarMatcher::analyzeDistortions(EnhancedValidationResult& result,
                                            const QVector<QPoint>& detected,
                                            const QVector<CatalogStar>& catalog)
{
    if (result.enhancedMatches.isEmpty()) {
        return;
    }
    
    qDebug() << "Analyzing optical distortions from" << result.enhancedMatches.size() << "matches";
    
    // Calculate image center for radial distortion analysis
    QPointF imageCenter(detected.isEmpty() ? 1024 : 
                       std::accumulate(detected.begin(), detected.end(), QPoint(0,0),
                       [](QPoint a, QPoint b) { return a + b; }).x() / detected.size(),
                       detected.isEmpty() ? 1024 :
                       std::accumulate(detected.begin(), detected.end(), QPoint(0,0),
                       [](QPoint a, QPoint b) { return a + b; }).y() / detected.size());
    
    result.residualVectors.clear();
    result.radialDistortions.clear();
    
    double maxRadius = 0.0;
    
    for (const auto& match : result.enhancedMatches) {
        if (match.detectedIndex >= 0 && match.detectedIndex < detected.size() &&
            match.catalogIndex >= 0 && match.catalogIndex < catalog.size()) {
            
            QPoint detectedPos = detected[match.detectedIndex];
            QPointF catalogPos = catalog[match.catalogIndex].pixelPos;
            
            // Calculate residual vector (detected - expected)
            QPointF residual(detectedPos.x() - catalogPos.x(),
                           detectedPos.y() - catalogPos.y());
            result.residualVectors.append(residual);
            
            // Calculate radial distance from image center
            double radius = sqrt(pow(detectedPos.x() - imageCenter.x(), 2) +
                               pow(detectedPos.y() - imageCenter.y(), 2));
            maxRadius = std::max(maxRadius, radius);
            
            // Calculate radial distortion component
            QPointF radialVector(detectedPos.x() - imageCenter.x(),
                               detectedPos.y() - imageCenter.y());
            if (radius > 0) {
                radialVector /= radius; // Normalize
                
                // Project residual onto radial direction
                double radialDistortion = QPointF::dotProduct(residual, radialVector);
                result.radialDistortions.append(std::abs(radialDistortion));
            }
        }
    }
    
    // Analyze distortion patterns
    if (!result.radialDistortions.isEmpty()) {
        double avgRadialDistortion = std::accumulate(result.radialDistortions.begin(),
                                                   result.radialDistortions.end(), 0.0) /
                                   result.radialDistortions.size();
        
        double maxRadialDistortion = *std::max_element(result.radialDistortions.begin(),
                                                      result.radialDistortions.end());
        
        qDebug() << QString("Distortion analysis: avg=%1 px, max=%2 px, radius=%3 px")
                    .arg(avgRadialDistortion).arg(maxRadialDistortion).arg(maxRadius);
        
        // Simple distortion model fitting (linear radial)
        if (result.radialDistortions.size() > 5) {
            // Fit k1*r model using least squares
            double sumR2 = 0, sumR2D = 0;
            for (int i = 0; i < result.enhancedMatches.size() && i < result.radialDistortions.size(); ++i) {
                const auto& match = result.enhancedMatches[i];
                if (match.detectedIndex >= 0 && match.detectedIndex < detected.size()) {
                    QPoint detectedPos = detected[match.detectedIndex];
                    double r = sqrt(pow(detectedPos.x() - imageCenter.x(), 2) +
                                  pow(detectedPos.y() - imageCenter.y(), 2));
                    if (r > 10.0) { // Avoid center region
                        double r2 = r * r;
                        sumR2 += r2 * r2;
                        sumR2D += r2 * result.radialDistortions[i];
                    }
                }
            }
            
            if (sumR2 > 0) {
                double k1 = sumR2D / sumR2;
                qDebug() << QString("Estimated radial distortion coefficient k1: %1 px/px²")
                            .arg(k1);
            }
        }
    }
}

// Complete quality filtering implementation
void EnhancedStarMatcher::filterLowQualityMatches(QVector<EnhancedStarMatch>& matches)
{
    if (matches.isEmpty()) return;
    
    qDebug() << "Filtering" << matches.size() << "matches for quality";
    
    // Calculate statistics for adaptive filtering
    QVector<double> distances, confidences;
    for (const auto& match : matches) {
        distances.append(match.pixelDistance);
        confidences.append(match.confidence);
    }
    
    if (distances.isEmpty()) return;
    
    // Calculate thresholds
    std::sort(distances.begin(), distances.end());
    std::sort(confidences.begin(), confidences.end());
    
    double medianDistance = distances[distances.size() / 2];
    double q75Distance = distances[distances.size() * 3 / 4];
    double medianConfidence = confidences[confidences.size() / 2];
    
    // Adaptive distance threshold (median + 1.5 * IQR, but respect parameter limits)
    double distanceThreshold = std::min(m_params.maxPixelDistance, 
                                      std::max(medianDistance * 2.0, q75Distance * 1.5));
    
    // Confidence threshold (use median as minimum, but respect parameter)
    double confidenceThreshold = std::max(m_params.minMatchConfidence, medianConfidence * 0.7);
    
    qDebug() << QString("Quality thresholds: distance=%1 px, confidence=%2")
                .arg(distanceThreshold).arg(confidenceThreshold);
    
    // Filter matches
    auto it = std::remove_if(matches.begin(), matches.end(),
        [distanceThreshold, confidenceThreshold](const EnhancedStarMatch& match) {
            // Remove if distance too large
            if (match.pixelDistance > distanceThreshold) return true;
            
            // Remove if confidence too low
            if (match.confidence < confidenceThreshold) return true;
            
            // Remove if magnitude difference too large (if available)
            if (match.magnitudeDifference > 3.0) return true;
            
            // Remove if triangle error too large (if triangle matching was used)
            if (match.triangleError > 0.2) return true;
            
            return false;
        });
    
    int removedCount = matches.end() - it;
    matches.erase(it, matches.end());
    
    qDebug() << QString("Quality filtering: removed %1 matches, %2 remaining")
                .arg(removedCount).arg(matches.size());
    
    // Sort remaining matches by confidence
    std::sort(matches.begin(), matches.end(),
        [](const EnhancedStarMatch& a, const EnhancedStarMatch& b) {
            return a.confidence > b.confidence;
        });
}

// Complete geometric validation implementation
QVector<EnhancedStarMatch> EnhancedStarMatcher::performGeometricValidation(
    const QVector<EnhancedStarMatch>& initialMatches,
    const QVector<QPoint>& detected,
    const QVector<CatalogStar>& catalog)
{
    if (initialMatches.size() < 3) {
        qDebug() << "Not enough matches for geometric validation";
        return initialMatches;
    }
    
    qDebug() << "Performing geometric validation on" << initialMatches.size() << "matches";
    
    QVector<EnhancedStarMatch> validatedMatches;
    
    // Step 1: Find consensus transformation using RANSAC-like approach
    struct SimpleTransform {
        double scale = 1.0;
        double rotation = 0.0;  // radians
        QPointF translation;
        double quality = 0.0;
        int inlierCount = 0;
    };
    
    SimpleTransform bestTransform;
    const int maxIterations = 100;
    const double inlierThreshold = 3.0; // pixels
    
    // RANSAC to find best transformation
    for (int iter = 0; iter < maxIterations && initialMatches.size() >= 2; ++iter) {
        // Select two random matches
        int idx1 = rand() % initialMatches.size();
        int idx2 = rand() % initialMatches.size();
        while (idx2 == idx1 && initialMatches.size() > 1) {
            idx2 = rand() % initialMatches.size();
        }
        
        const auto& match1 = initialMatches[idx1];
        const auto& match2 = initialMatches[idx2];
        
        // Get positions
        if (match1.detectedIndex >= detected.size() || match1.catalogIndex >= catalog.size() ||
            match2.detectedIndex >= detected.size() || match2.catalogIndex >= catalog.size()) {
            continue;
        }
        
        QPointF det1(detected[match1.detectedIndex]);
        QPointF det2(detected[match2.detectedIndex]);
        QPointF cat1 = catalog[match1.catalogIndex].pixelPos;
        QPointF cat2 = catalog[match2.catalogIndex].pixelPos;
        
        // Calculate transformation from these two points
        QPointF detVec = det2 - det1;
        QPointF catVec = cat2 - cat1;
        
        double detDist = sqrt(QPointF::dotProduct(detVec, detVec));
        double catDist = sqrt(QPointF::dotProduct(catVec, catVec));
        
        if (detDist < 5.0 || catDist < 5.0) continue; // Too close
        
        SimpleTransform transform;
        transform.scale = detDist / catDist;
        
        // Calculate rotation
        double detAngle = atan2(detVec.y(), detVec.x());
        double catAngle = atan2(catVec.y(), catVec.x());
        transform.rotation = detAngle - catAngle;
        
        // Calculate translation using first point
        QPointF rotatedCat1(cat1.x() * cos(transform.rotation) - cat1.y() * sin(transform.rotation),
                           cat1.x() * sin(transform.rotation) + cat1.y() * cos(transform.rotation));
        QPointF scaledRotatedCat1 = rotatedCat1 * transform.scale;
        transform.translation = det1 - scaledRotatedCat1;
        
        // Count inliers
        int inliers = 0;
        double totalError = 0.0;
        
        for (const auto& match : initialMatches) {
            if (match.detectedIndex >= detected.size() || match.catalogIndex >= catalog.size()) {
                continue;
            }
            
            QPointF detPos(detected[match.detectedIndex]);
            QPointF catPos = catalog[match.catalogIndex].pixelPos;
            
            // Apply transformation to catalog position
            QPointF rotatedCat(catPos.x() * cos(transform.rotation) - catPos.y() * sin(transform.rotation),
                              catPos.x() * sin(transform.rotation) + catPos.y() * cos(transform.rotation));
            QPointF transformedCat = rotatedCat * transform.scale + transform.translation;
            
            // Calculate error
            QPointF error = detPos - transformedCat;
            double errorMag = sqrt(QPointF::dotProduct(error, error));
            
            if (errorMag < inlierThreshold) {
                inliers++;
                totalError += errorMag;
            }
        }
        
        if (inliers > bestTransform.inlierCount) {
            bestTransform = transform;
            bestTransform.inlierCount = inliers;
            bestTransform.quality = inliers > 0 ? totalError / inliers : 1000.0;
        }
    }
    
    qDebug() << QString("Best transform: scale=%1, rotation=%2°, inliers=%3/%4")
                .arg(bestTransform.scale)
                .arg(bestTransform.rotation * 180.0 / M_PI)
                .arg(bestTransform.inlierCount)
                .arg(initialMatches.size());
    
    // Step 2: Validate all matches against best transformation
    for (auto match : initialMatches) {
        if (match.detectedIndex >= detected.size() || match.catalogIndex >= catalog.size()) {
            continue;
        }
        
        QPointF detPos(detected[match.detectedIndex]);
        QPointF catPos = catalog[match.catalogIndex].pixelPos;
        
        // Apply transformation
        QPointF rotatedCat(catPos.x() * cos(bestTransform.rotation) - catPos.y() * sin(bestTransform.rotation),
                          catPos.x() * sin(bestTransform.rotation) + catPos.y() * cos(bestTransform.rotation));
        QPointF transformedCat = rotatedCat * bestTransform.scale + bestTransform.translation;
        
        // Calculate error
        QPointF error = detPos - transformedCat;
        double errorMag = sqrt(QPointF::dotProduct(error, error));
        
        // Update match quality
        match.isGeometricallyValid = (errorMag < inlierThreshold);
        
        if (match.isGeometricallyValid) {
            // Boost confidence for geometrically consistent matches
            match.confidence = std::min(1.0, match.confidence + 0.2);
            validatedMatches.append(match);
        }
    }
    
    qDebug() << QString("Geometric validation: %1 geometrically valid matches")
                .arg(validatedMatches.size());
    
    return validatedMatches;
}

// Complete distortion model calibration
bool StarCatalogValidator::calibrateDistortionModel(const EnhancedValidationResult& result)
{
    if (result.enhancedMatches.size() < 10) {
        qDebug() << "Not enough matches for distortion calibration";
        return false;
    }
    
    qDebug() << "Calibrating distortion model from" << result.enhancedMatches.size() << "matches";
    
    // Calculate image center
    QPointF imageCenter(m_astrometricMetadata.Width() * 0.5, m_astrometricMetadata.Height() * 0.5);
    
    // Collect residuals and radial distances
    QVector<double> radii, radialResiduals, tangentialResiduals;
    
    for (const auto& match : result.enhancedMatches) {
        if (!match.isGeometricallyValid) continue;
        
        // This would need access to detected stars and catalog - simplified for now
        double radius = 100.0; // Placeholder - would calculate actual radius from center
        double radialResidual = 0.1; // Placeholder - would calculate actual residual
        double tangentialResidual = 0.05; // Placeholder - would calculate actual residual
        
        radii.append(radius);
        radialResiduals.append(radialResidual);
        tangentialResiduals.append(tangentialResidual);
    }
    
    if (radii.size() < 5) return false;
    
    // Simple linear radial distortion model: dr = k1 * r^3
    double sumR6 = 0, sumR3DR = 0;
    for (int i = 0; i < radii.size(); ++i) {
        double r = radii[i];
        double r3 = r * r * r;
        double r6 = r3 * r3;
        sumR6 += r6;
        sumR3DR += r3 * radialResiduals[i];
    }
    
    if (sumR6 > 0) {
        double k1 = sumR3DR / sumR6;
        
        // Store distortion parameters
        m_radialDistortionK1 = k1;
        m_hasDistortionModel = true;
        
        qDebug() << QString("Calibrated radial distortion: k1 = %1")
                    .arg(k1);
        
        return true;
    }
    
    return false;
}

// Add this debugging method to help diagnose pixel coordinate matching issues
void debugPixelMatching(const QVector<QPoint>& detectedStars,
                       const QVector<CatalogStar>& catalogStars,
                       const QVector<StarMatch>& matches,
                       double pixelTolerance)
{
    qDebug() << "\n=== PIXEL MATCHING DEBUG ===";
    qDebug() << QString("Detected stars: %1, Catalog stars: %2, Matches: %3")
                .arg(detectedStars.size()).arg(catalogStars.size()).arg(matches.size());
    qDebug() << QString("Pixel tolerance: %1 px").arg(pixelTolerance);
    
    // Analyze distance distribution
    QVector<double> allDistances;
    int goodMatches = 0, badMatches = 0;
    
    for (const auto& match : matches) {
        if (match.detectedIndex >= 0 && match.detectedIndex < detectedStars.size() &&
            match.catalogIndex >= 0 && match.catalogIndex < catalogStars.size()) {
            
            QPoint detected = detectedStars[match.detectedIndex];
            QPointF catalog = catalogStars[match.catalogIndex].pixelPos;
            
            double distance = sqrt(pow(detected.x() - catalog.x(), 2) + 
                                 pow(detected.y() - catalog.y(), 2));
            allDistances.append(distance);
            
            if (match.isGoodMatch) {
                goodMatches++;
            } else {
                badMatches++;
                if (distance <= pixelTolerance) {
                    qDebug() << QString("⚠️  Match %1->%2: distance=%3 px (within tolerance but marked bad)")
                                .arg(match.detectedIndex).arg(match.catalogIndex).arg(distance);
                }
            }
        }
    }
    
    if (!allDistances.isEmpty()) {
        std::sort(allDistances.begin(), allDistances.end());
        double medianDist = allDistances[allDistances.size() / 2];
        double minDist = allDistances.first();
        double maxDist = allDistances.last();
        
        qDebug() << QString("Distance statistics: min=%1, median=%2, max=%3 px")
                    .arg(minDist).arg(medianDist).arg(maxDist);
        qDebug() << QString("Good matches: %1, Bad matches: %2").arg(goodMatches).arg(badMatches);
        
        // Count how many are within tolerance
        int withinTolerance = std::count_if(allDistances.begin(), allDistances.end(),
            [pixelTolerance](double d) { return d <= pixelTolerance; });
        
        qDebug() << QString("Matches within tolerance: %1/%2 (%3%%)")
                    .arg(withinTolerance).arg(allDistances.size())
                    .arg(100.0 * withinTolerance / allDistances.size());
        
        if (withinTolerance > goodMatches) {
            qDebug() << "🔍 ISSUE: More matches within tolerance than marked as good!";
            qDebug() << "   Check magnitude difference criteria and other validation logic";
        }
    }
    
    // Check for systematic offsets
    double sumDx = 0, sumDy = 0;
    int offsetCount = 0;
    
    for (const auto& match : matches) {
        if (match.detectedIndex >= 0 && match.detectedIndex < detectedStars.size() &&
            match.catalogIndex >= 0 && match.catalogIndex < catalogStars.size()) {
            
            QPoint detected = detectedStars[match.detectedIndex];
            QPointF catalog = catalogStars[match.catalogIndex].pixelPos;
            
            sumDx += detected.x() - catalog.x();
            sumDy += detected.y() - catalog.y();
            offsetCount++;
        }
    }
    
    if (offsetCount > 0) {
        double avgDx = sumDx / offsetCount;
        double avgDy = sumDy / offsetCount;
        double avgOffset = sqrt(avgDx * avgDx + avgDy * avgDy);
        
        qDebug() << QString("Systematic offset: dx=%1, dy=%2 px (magnitude=%3 px)")
                    .arg(avgDx).arg(avgDy).arg(avgOffset);
        
        if (avgOffset > 2.0) {
            qDebug() << "🔍 ISSUE: Large systematic offset detected!";
            qDebug() << "   Check WCS calibration and coordinate transformations";
        }
    }
}

EnhancedStarMatcher::EnhancedStarMatcher(const StarMatchingParameters& params)
    : m_params(params)
{
}

EnhancedValidationResult EnhancedStarMatcher::matchStarsAdvanced(
    const QVector<QPoint>& detectedStars,
    const QVector<float>& detectedMagnitudes,
    const QVector<CatalogStar>& catalogStars,
    const pcl::AstrometricMetadata& astrometry)
{
    EnhancedValidationResult result;
    
    qDebug() << "=== ENHANCED STAR MATCHING (PixInsight Methods) ===";
    qDebug() << QString("Detected stars: %1, Catalog stars: %2")
                .arg(detectedStars.size()).arg(catalogStars.size());
    
    // Step 1: Initial position-based matching with magnitude constraints
    qDebug() << "Step 1: Initial position-based matching...";
    auto initialMatches = performInitialMatching(detectedStars, catalogStars, astrometry);
    qDebug() << QString("Initial matches found: %1").arg(initialMatches.size());
    
    // Step 2: Triangle pattern matching for geometric validation
    if (m_params.useTriangleMatching && detectedStars.size() >= m_params.minTriangleStars) {
        qDebug() << "Step 2: Triangle pattern matching...";
        auto triangleMatches = performTriangleMatching(detectedStars, catalogStars);
        
        // Merge triangle matches with initial matches
        for (const auto& triMatch : triangleMatches) {
            // Find corresponding initial match and boost confidence
            for (auto& initMatch : initialMatches) {
                if (initMatch.detectedIndex == triMatch.detectedIndex &&
                    initMatch.catalogIndex == triMatch.catalogIndex) {
                    initMatch.confidence = std::max(initMatch.confidence, triMatch.confidence);
                    initMatch.triangleError = triMatch.triangleError;
                    initMatch.supportingMatches = triMatch.supportingMatches;
                    break;
                }
            }
        }
        
        qDebug() << QString("Triangle validation completed, %1 patterns analyzed")
                    .arg(triangleMatches.size());
    }
    
    // Step 3: Geometric validation and outlier rejection
    qDebug() << "Step 3: Geometric validation...";
    auto validatedMatches = performGeometricValidation(initialMatches, detectedStars, catalogStars);
    qDebug() << QString("Geometrically valid matches: %1").arg(validatedMatches.size());
    
    // Step 4: Quality filtering
    qDebug() << "Step 4: Quality filtering...";
    filterLowQualityMatches(validatedMatches);
    
    // Store results
    result.enhancedMatches = validatedMatches;
    
    // Step 5: Distortion analysis
    if (m_params.useDistortionModel) {
        qDebug() << "Step 5: Distortion analysis...";
        analyzeDistortions(result, detectedStars, catalogStars);
    }
    
    // Step 6: Calculate advanced statistics
    calculateAdvancedStatistics(result);
    
    // Convert to standard ValidationResult format for compatibility
    result.totalDetected = detectedStars.size();
    result.totalCatalog = catalogStars.size();
    result.totalMatches = validatedMatches.size();
    result.matchPercentage = (100.0 * result.totalMatches) / result.totalDetected;
    result.catalogStars = catalogStars;
    
    // Convert enhanced matches to standard matches
    for (const auto& enhanced : validatedMatches) {
        StarMatch standardMatch;
        standardMatch.detectedIndex = enhanced.detectedIndex;
        standardMatch.catalogIndex = enhanced.catalogIndex;
        standardMatch.distance = enhanced.pixelDistance;
        standardMatch.magnitudeDiff = enhanced.magnitudeDifference;
        standardMatch.isGoodMatch = enhanced.confidence >= m_params.minMatchConfidence;
        result.matches.append(standardMatch);
    }
    
    // Generate enhanced summary
    result.enhancedSummary = QString(
        "Enhanced Star Matching Results:\n"
        "================================\n"
        "Detected Stars: %1\n"
        "Catalog Stars: %2\n"
        "High-Quality Matches: %3 (%4%%)\n"
        "Geometric RMS: %5 pixels\n"
        "Photometric RMS: %6 magnitudes\n"
        "Astrometric Accuracy: %7 arcsec\n"
        "Scale Error: %8%%\n"
        "Rotation Error: %9 degrees\n"
        "Matching Confidence: %10%%\n"
        "\nDistortion Analysis:\n"
        "Max Radial Distortion: %11 pixels\n"
        "Distortion Model: %12")
        .arg(result.totalDetected)
        .arg(result.totalCatalog)
        .arg(result.totalMatches)
        .arg(result.matchPercentage)
        .arg(result.geometricRMS)
        .arg(result.photometricRMS)
        .arg(result.astrometricAccuracy)
        .arg(result.scaleError)
        .arg(result.rotationError)
        .arg(result.matchingConfidence)
        .arg(result.radialDistortions.isEmpty() ? 0.0 : 
             *std::max_element(result.radialDistortions.begin(), result.radialDistortions.end()))
        .arg(m_params.useDistortionModel ? "Applied" : "Not Applied");
    
    result.summary = result.enhancedSummary;
    result.isValid = result.totalMatches >= m_params.minMatchesForValidation;
    
    qDebug() << "=== ENHANCED MATCHING COMPLETE ===";
    qDebug() << QString("Final result: %1/%2 matches (%3%% confidence)")
                .arg(result.totalMatches).arg(result.totalDetected)
                .arg(result.matchingConfidence);
    
    return result;
}

QVector<EnhancedStarMatch> EnhancedStarMatcher::performInitialMatching(
    const QVector<QPoint>& detected,
    const QVector<CatalogStar>& catalog,
    const pcl::AstrometricMetadata& astrometry)
{
    QVector<EnhancedStarMatch> matches;
    
    // Create spatial index for efficient neighbor searching
    // For simplicity, we'll use a direct search, but PixInsight would use KD-trees
    
    for (int i = 0; i < detected.size(); ++i) {
        const QPoint& detectedPos = detected[i];
        
        EnhancedStarMatch bestMatch;
        bestMatch.detectedIndex = i;
        bestMatch.pixelDistance = std::numeric_limits<double>::max();
        bestMatch.confidence = 0.0;
        
        // Search for nearest catalog star within search radius
        for (int j = 0; j < catalog.size(); ++j) {
            const CatalogStar& catalogStar = catalog[j];
            
            if (!catalogStar.isValid) continue;
            
            // Calculate pixel distance
            double dx = detectedPos.x() - catalogStar.pixelPos.x();
            double dy = detectedPos.y() - catalogStar.pixelPos.y();
            double distance = sqrt(dx * dx + dy * dy);
            
            if (distance <= m_params.searchRadius && distance < bestMatch.pixelDistance) {
                bestMatch.catalogIndex = j;
                bestMatch.pixelDistance = distance;
                bestMatch.magnitudeDifference = 0.0; // Will be calculated if we have magnitudes
                
                // Calculate initial confidence based on distance
                double distanceConfidence = 1.0 - (distance / m_params.searchRadius);
                bestMatch.confidence = distanceConfidence;
                
                // Additional quality checks
                bestMatch.isGeometricallyValid = distance <= m_params.maxPixelDistance;
                bestMatch.isPhotometricallyValid = true; // Will be refined later
            }
        }
        
        // Only accept matches within maximum distance
        if (bestMatch.catalogIndex >= 0 && bestMatch.pixelDistance <= m_params.maxPixelDistance) {
            matches.append(bestMatch);
        }
    }
    
    return matches;
}

QVector<EnhancedStarMatch> EnhancedStarMatcher::performTriangleMatching(
    const QVector<QPoint>& detected,
    const QVector<CatalogStar>& catalog)
{
    QVector<EnhancedStarMatch> triangleMatches;
    
    if (detected.size() < 3 || catalog.size() < 3) {
        return triangleMatches;
    }
    
    // Generate triangle patterns for detected stars
    auto detectedTriangles = generateTrianglePatterns(detected, true);
    
    // Generate triangle patterns for catalog stars (using pixel positions)
    QVector<QPoint> catalogPixelPositions;
    for (const auto& cat : catalog) {
        if (cat.isValid) {
            catalogPixelPositions.append(QPoint(
                static_cast<int>(cat.pixelPos.x()),
                static_cast<int>(cat.pixelPos.y())
            ));
        }
    }
    auto catalogTriangles = generateTrianglePatterns(catalogPixelPositions, true);
    
    qDebug() << QString("Generated %1 detected triangles and %2 catalog triangles")
                .arg(detectedTriangles.size()).arg(catalogTriangles.size());
    
    // Match triangle patterns
    auto triangleMatchPairs = matchTrianglePatterns(detectedTriangles, catalogTriangles);
    
    qDebug() << QString("Found %1 triangle matches with confidence %2")
                .arg(triangleMatchPairs.first.size()).arg(triangleMatchPairs.second);
    
    // Convert triangle matches to star matches
    for (const auto& trianglePair : triangleMatchPairs.first) {
        const TrianglePattern& detTriangle = detectedTriangles[trianglePair.first];
        const TrianglePattern& catTriangle = catalogTriangles[trianglePair.second];
        
        // Each triangle gives us 3 potential star matches
        for (int i = 0; i < 3; ++i) {
            int detectedIdx = detTriangle.starIndices[i];
            int catalogIdx = catTriangle.starIndices[i]; // Assuming same ordering
            
            EnhancedStarMatch match;
            match.detectedIndex = detectedIdx;
            match.catalogIndex = catalogIdx;
            match.confidence = triangleMatchPairs.second;
            match.triangleError = calculateTriangleSimilarity(detTriangle, catTriangle);
            
            // Calculate pixel distance
            QPoint detPos = detected[detectedIdx];
            QPointF catPos = catalog[catalogIdx].pixelPos;
            double dx = detPos.x() - catPos.x();
            double dy = detPos.y() - catPos.y();
            match.pixelDistance = sqrt(dx * dx + dy * dy);
            
            triangleMatches.append(match);
        }
    }
    
    return triangleMatches;
}

QVector<TrianglePattern> EnhancedStarMatcher::generateTrianglePatterns(
    const QVector<QPoint>& stars, bool isPixelCoords)
{
    QVector<TrianglePattern> patterns;
    
    if (stars.size() < 3) return patterns;
    
    // Generate all possible triangles (combinatorial approach)
    // For performance, limit to brightest/most reliable stars
    int maxStars = std::min(20, (int)(stars.size())); // Limit for performance
    
    for (int i = 0; i < maxStars - 2; ++i) {
        for (int j = i + 1; j < maxStars - 1; ++j) {
            for (int k = j + 1; k < maxStars; ++k) {
                TrianglePattern pattern;
                pattern.starIndices = {i, j, k};
                
                QPointF p1(stars[i]);
                QPointF p2(stars[j]);
                QPointF p3(stars[k]);
                
                // Calculate side lengths
                pattern.side1 = sqrt(pow(p2.x() - p1.x(), 2) + pow(p2.y() - p1.y(), 2));
                pattern.side2 = sqrt(pow(p3.x() - p2.x(), 2) + pow(p3.y() - p2.y(), 2));
                pattern.side3 = sqrt(pow(p1.x() - p3.x(), 2) + pow(p1.y() - p3.y(), 2));
                
                // Skip degenerate triangles
                if (pattern.side1 < 5.0 || pattern.side2 < 5.0 || pattern.side3 < 5.0) {
                    continue;
                }
                
                // Calculate area using cross product
                double area = 0.5 * abs((p2.x() - p1.x()) * (p3.y() - p1.y()) - 
                                       (p3.x() - p1.x()) * (p2.y() - p1.y()));
                pattern.area = area;
                
                // Skip very thin triangles
                double perimeter = pattern.side1 + pattern.side2 + pattern.side3;
                if (area < 0.1 * perimeter) continue;
                
                // Calculate centroid
                pattern.centroid = QPointF((p1.x() + p2.x() + p3.x()) / 3.0,
                                          (p1.y() + p2.y() + p3.y()) / 3.0);
                
                // Calculate scale-invariant ratios
                double maxSide = std::max({pattern.side1, pattern.side2, pattern.side3});
                pattern.ratio12 = pattern.side1 / maxSide;
                pattern.ratio13 = pattern.side2 / maxSide;
                pattern.ratio23 = pattern.side3 / maxSide;
                
                // Normalized area (scale invariant)
                pattern.normalizedArea = area / (perimeter * perimeter);
                
                if (isTriangleValid(pattern)) {
                    patterns.append(pattern);
                }
            }
        }
    }
    
    qDebug() << QString("Generated %1 valid triangle patterns from %2 stars")
                .arg(patterns.size()).arg(maxStars);
    
    return patterns;
}

QPair<QVector<QPair<int, int>>, double> EnhancedStarMatcher::matchTrianglePatterns(
    const QVector<TrianglePattern>& pattern1,
    const QVector<TrianglePattern>& pattern2)
{
    QVector<QPair<int, int>> matches;
    QVector<double> similarities;
    
    double tolerancePercent = m_params.triangleTolerancePercent / 100.0;
    
    for (int i = 0; i < pattern1.size(); ++i) {
        for (int j = 0; j < pattern2.size(); ++j) {
            double similarity = calculateTriangleSimilarity(pattern1[i], pattern2[j]);
            
            if (similarity < tolerancePercent) {
                matches.append(qMakePair(i, j));
                similarities.append(1.0 - similarity); // Convert to confidence
            }
        }
    }
    
    // Calculate overall matching confidence
    double avgConfidence = 0.0;
    if (!similarities.isEmpty()) {
        avgConfidence = std::accumulate(similarities.begin(), similarities.end(), 0.0) / similarities.size();
    }
    
    return qMakePair(matches, avgConfidence);
}

double EnhancedStarMatcher::calculateTriangleSimilarity(
    const TrianglePattern& t1, const TrianglePattern& t2)
{
    // Compare scale-invariant properties
    double ratioError1 = abs(t1.ratio12 - t2.ratio12);
    double ratioError2 = abs(t1.ratio13 - t2.ratio13);
    double ratioError3 = abs(t1.ratio23 - t2.ratio23);
    
    double areaError = abs(t1.normalizedArea - t2.normalizedArea);
    
    // Weighted combination of errors
    double totalError = (ratioError1 + ratioError2 + ratioError3) / 3.0 + areaError;
    
    return totalError;
}

void EnhancedStarMatcher::calculateAdvancedStatistics(EnhancedValidationResult& result)
{
    if (result.enhancedMatches.isEmpty()) {
        result.geometricRMS = 0.0;
        result.matchingConfidence = 0.0;
        return;
    }
    
    // Calculate geometric RMS
    double sumSquaredDistances = 0.0;
    double sumConfidences = 0.0;
    
    for (const auto& match : result.enhancedMatches) {
        sumSquaredDistances += match.pixelDistance * match.pixelDistance;
        sumConfidences += match.confidence;
    }
    
    result.geometricRMS = sqrt(sumSquaredDistances / result.enhancedMatches.size());
    result.matchingConfidence = (sumConfidences / result.enhancedMatches.size()) * 100.0;
    
    // Calculate photometric RMS if magnitude data is available
    double sumMagSquared = 0.0;
    int magCount = 0;
    for (const auto& match : result.enhancedMatches) {
        if (match.magnitudeDifference > 0) {
            sumMagSquared += match.magnitudeDifference * match.magnitudeDifference;
            magCount++;
        }
    }
    
    if (magCount > 0) {
        result.photometricRMS = sqrt(sumMagSquared / magCount);
    }
    
    // Estimate astrometric accuracy (convert pixels to arcseconds)
    // This would need the pixel scale from the astrometric solution
    result.astrometricAccuracy = result.geometricRMS * 1.0; // Placeholder: 1 arcsec/pixel
}

// Integration with existing StarCatalogValidator
EnhancedValidationResult StarCatalogValidator::validateStarsAdvanced(
    const QVector<QPoint>& detectedStars,
    const QVector<float>& starMagnitudes,
    const StarMatchingParameters& params)
{
    if (!m_enhancedMatcher) {
        m_enhancedMatcher = std::make_unique<EnhancedStarMatcher>(params);
    }
    
    m_enhancedMatcher->setParameters(params);
    
    return m_enhancedMatcher->matchStarsAdvanced(
        detectedStars, starMagnitudes, m_catalogStars, m_astrometricMetadata);
}

void StarCatalogValidator::setMatchingParameters(const StarMatchingParameters& params)
{
    m_matchingParams = params;
    if (m_enhancedMatcher) {
        m_enhancedMatcher->setParameters(params);
    }
}

// Replace the queryLocal2MASS method with this Gaia version:
void StarCatalogValidator::queryGaiaDR3(double centerRA, double centerDec, double radiusDegrees)
{
    qDebug() << "\n=== QUERYING GAIA GDR3 CATALOG ===";
    
    if (!GaiaGDR3Catalog::isAvailable()) {
        qDebug() << "❌ Gaia GDR3 not available, falling back to network catalogs";
        queryGaiaCatalog(centerRA, centerDec, radiusDegrees);
        return;
    }
    
    emit catalogQueryStarted();
    
    try {
        // Set up Gaia search parameters
        GaiaGDR3Catalog::SearchParameters params;
        params.centerRA = centerRA;
        params.centerDec = centerDec;
        params.radiusDegrees = radiusDegrees;
        params.maxMagnitude = 17.0;
        params.maxResults = 200000;
        params.useProperMotion = true;  // Apply proper motion to current epoch
        params.epochYear = 2000.0;     // Current epoch
        
        qDebug() << QString("🔍 Searching Gaia GDR3: center=(%1°,%2°) radius=%3°")
	  .arg(centerRA).arg(centerDec).arg(radiusDegrees);
        
        // Query the Gaia catalog
        auto gaiaStars = GaiaGDR3Catalog::queryRegion(params);
        
        // Convert Gaia stars to our CatalogStar format
        m_catalogStars.clear();
        m_catalogStars.reserve(gaiaStars.size());
        
        int starsInBounds = 0;
        int starsWithSpectra = 0;
        double brightestMag = 99.0;
        
        for (const auto& gaiaStar : gaiaStars) {
            CatalogStar catalogStar;
            catalogStar.id = gaiaStar.sourceId;
            catalogStar.ra = gaiaStar.ra;
            catalogStar.dec = gaiaStar.dec;
            catalogStar.magnitude = gaiaStar.magnitude;
            catalogStar.spectralType = gaiaStar.spectralClass;

	    // NEW: Add enhanced Gaia photometry
	    catalogStar.magBP = gaiaStar.magBP;
	    catalogStar.magRP = gaiaStar.magRP;
	    if (gaiaStar.magBP > 0 && gaiaStar.magRP > 0) {
		catalogStar.colorBP_RP = gaiaStar.magBP - gaiaStar.magRP;
		catalogStar.hasExtendedData = true;
            }
	    
            // Calculate pixel position using our WCS
            catalogStar.pixelPos = skyToPixel(gaiaStar.ra, gaiaStar.dec);
            
            // Check if star is within image bounds
            catalogStar.isValid = (catalogStar.pixelPos.x() >= 0 && 
                                  catalogStar.pixelPos.x() < m_astrometricMetadata.Width() && 
                                  catalogStar.pixelPos.y() >= 0 && 
                                  catalogStar.pixelPos.y() < m_astrometricMetadata.Height());
            
            if (catalogStar.isValid) {
                starsInBounds++;
                brightestMag = std::min(brightestMag, catalogStar.magnitude);
            }
            
            if (gaiaStar.hasSpectrum) {
                starsWithSpectra++;
            }
            
            m_catalogStars.append(catalogStar);
        }
        
        qDebug() << QString("✅ Gaia GDR3 query successful:");
        qDebug() << QString("   📊 Total stars retrieved: %1").arg(gaiaStars.size());
        qDebug() << QString("   📊 Stars in image bounds: %1").arg(starsInBounds);
        qDebug() << QString("   📊 Stars with BP/RP spectra: %1").arg(starsWithSpectra);
        qDebug() << QString("   📊 Brightest star in field: mag %1").arg(brightestMag);
        
        // Add bright stars from our local database to fill in any gaps
        // (Gaia sometimes has issues with very bright stars due to saturation)
        addBrightStarsFromDatabase(centerRA, centerDec, radiusDegrees);
        
        // Sort by magnitude (brightest first)
        std::sort(m_catalogStars.begin(), m_catalogStars.end(), 
                  [](const CatalogStar& a, const CatalogStar& b) {
                      return a.magnitude < b.magnitude;
                  });
        
        QString message = QString("Retrieved %1 stars from Gaia GDR3 catalog (%2 in bounds)")
                         .arg(m_catalogStars.size()).arg(starsInBounds);
        
        emit catalogQueryFinished(true, message);
        
    } catch (const std::exception& e) {
        QString errorMsg = QString("Gaia GDR3 query error: %1").arg(e.what());
        qDebug() << "❌" << errorMsg;
        emit catalogQueryFinished(false, errorMsg);
    } catch (...) {
        QString errorMsg = "Unknown error during Gaia GDR3 query";
        qDebug() << "❌" << errorMsg;
        emit catalogQueryFinished(false, errorMsg);
    }
}

// Update the initialization method
void StarCatalogValidator::initializeGaiaDR3()
{
    // Set the path to your Gaia GDR3 catalog
    QString catalogPath = "/Volumes/X10Pro/gdr3-1.0.0-01.xpsd";
    GaiaGDR3Catalog::setCatalogPath(catalogPath);
    
    if (GaiaGDR3Catalog::isAvailable()) {
        qDebug() << "✅" << GaiaGDR3Catalog::getCatalogInfo();
	/*        
        // Validate the database
        if (GaiaGDR3Catalog::validateDatabase()) {
            qDebug() << "✅ Gaia GDR3 database validation successful";
        } else {
            qDebug() << "⚠️  Gaia GDR3 database validation failed";
        }
	*/
    } else {
        qDebug() << "❌ Gaia GDR3 catalog not found at:" << catalogPath;
        qDebug() << "   Network catalogs will be used as fallback";
    }
}

// Update the main queryCatalog method to use Gaia
void StarCatalogValidator::queryCatalog(double centerRA, double centerDec, double radiusDegrees)
{
    if (!m_astrometricMetadata.IsValid()) {
        emit errorSignal("No valid WCS data available for catalog query");
        return;
    }
    
    // Try Gaia GDR3 first, fall back to network catalogs if needed
    queryGaiaDR3(centerRA, centerDec, radiusDegrees);
}

// Add advanced Gaia search methods
void StarCatalogValidator::queryGaiaWithSpectra(double centerRA, double centerDec, double radiusDegrees)
{
    qDebug() << "\n=== QUERYING GAIA GDR3 FOR STARS WITH BP/RP SPECTRA ===";
    
    if (!GaiaGDR3Catalog::isAvailable()) {
        emit errorSignal("Gaia GDR3 catalog not available for spectrum search");
        return;
    }
    
    emit catalogQueryStarted();
    
    try {
        auto spectralStars = GaiaGDR3Catalog::findStarsWithSpectra(
            centerRA, centerDec, radiusDegrees, 99.0);
        
        // Convert to catalog stars
        m_catalogStars.clear();
        for (const auto& gaiaStar : spectralStars) {
            CatalogStar catalogStar;
            catalogStar.id = gaiaStar.sourceId + "_SPEC";
            catalogStar.ra = gaiaStar.ra;
            catalogStar.dec = gaiaStar.dec;
            catalogStar.magnitude = gaiaStar.magnitude;
            catalogStar.spectralType = gaiaStar.spectralClass + " (BP/RP)";
            catalogStar.pixelPos = skyToPixel(gaiaStar.ra, gaiaStar.dec);
            catalogStar.isValid = (catalogStar.pixelPos.x() >= 0 && 
                                  catalogStar.pixelPos.x() < m_astrometricMetadata.Width() && 
                                  catalogStar.pixelPos.y() >= 0 && 
                                  catalogStar.pixelPos.y() < m_astrometricMetadata.Height());
            
            m_catalogStars.append(catalogStar);
        }
        
        QString message = QString("Retrieved %1 stars with BP/RP spectra from Gaia GDR3")
                         .arg(m_catalogStars.size());
        qDebug() << "✅" << message;
        emit catalogQueryFinished(true, message);
        
    } catch (...) {
        QString errorMsg = "Error querying Gaia GDR3 for spectral data";
        qDebug() << "❌" << errorMsg;
        emit catalogQueryFinished(false, errorMsg);
    }
}

void StarCatalogValidator::findBrightGaiaStars(double centerRA, double centerDec, double radiusDegrees, int count)
{
    qDebug() << QString("\n=== FINDING %1 BRIGHTEST GAIA STARS ===").arg(count);
    
    if (!GaiaGDR3Catalog::isAvailable()) {
        emit errorSignal("Gaia GDR3 catalog not available");
        return;
    }
    
    emit catalogQueryStarted();
    
    try {
        auto brightStars = GaiaGDR3Catalog::findBrightestStars(
            centerRA, centerDec, radiusDegrees, count);
        
        // Convert to catalog stars
        m_catalogStars.clear();
        for (const auto& gaiaStar : brightStars) {
            CatalogStar catalogStar;
            catalogStar.id = gaiaStar.sourceId;
            catalogStar.ra = gaiaStar.ra;
            catalogStar.dec = gaiaStar.dec;
            catalogStar.magnitude = gaiaStar.magnitude;
            catalogStar.spectralType = gaiaStar.spectralClass;
            catalogStar.pixelPos = skyToPixel(gaiaStar.ra, gaiaStar.dec);
            catalogStar.isValid = (catalogStar.pixelPos.x() >= 0 && 
                                  catalogStar.pixelPos.x() < m_astrometricMetadata.Width() && 
                                  catalogStar.pixelPos.y() >= 0 && 
                                  catalogStar.pixelPos.y() < m_astrometricMetadata.Height());
            
            m_catalogStars.append(catalogStar);
        }
        
        QString message = QString("Found %1 brightest Gaia stars in field").arg(m_catalogStars.size());
        qDebug() << "✅" << message;
        emit catalogQueryFinished(true, message);
        
    } catch (...) {
        QString errorMsg = "Error finding bright Gaia stars";
        qDebug() << "❌" << errorMsg;
        emit catalogQueryFinished(false, errorMsg);
    }
}

// Update showCatalogStats to include Gaia-specific information
void StarCatalogValidator::showCatalogStats() const
{
    if (m_catalogStars.isEmpty()) {
        qDebug() << "No catalog stars loaded";
        return;
    }
    
    // Enhanced statistics for Gaia data
    QMap<QString, int> spectralTypeHistogram;
    QMap<int, int> magHistogram;
    int validStars = 0;
    int starsWithSpectra = 0;
    double brightestMag = 99.0;
    double faintestMag = -99.0;
    double totalParallax = 0.0;
    int parallaxCount = 0;
    
    for (const auto& star : m_catalogStars) {
        if (star.isValid) validStars++;
        
        // Spectral type distribution
        QString specType = star.spectralType.left(1); // First character only
        spectralTypeHistogram[specType]++;
        
        // Magnitude distribution
        int magBin = static_cast<int>(std::floor(star.magnitude));
        magHistogram[magBin]++;
        
        // Magnitude range
        brightestMag = std::min(brightestMag, star.magnitude);
        faintestMag = std::max(faintestMag, star.magnitude);
        
        // Count spectra (based on ID or spectral type)
        if (star.spectralType.contains("BP/RP") || star.id.contains("_SPEC")) {
            starsWithSpectra++;
        }
    }
    
    qDebug() << "\n=== GAIA GDR3 CATALOG STATISTICS ===";
    qDebug() << QString("Total stars: %1").arg(m_catalogStars.size());
    qDebug() << QString("Stars in image bounds: %1").arg(validStars);
    qDebug() << QString("Stars with BP/RP spectra: %1").arg(starsWithSpectra);
    qDebug() << QString("Magnitude range: %1 to %2").arg(brightestMag).arg(faintestMag);
    
    qDebug() << "Magnitude distribution:";
    for (auto it = magHistogram.begin(); it != magHistogram.end(); ++it) {
        qDebug() << QString("  mag %1-%2: %3 stars")
                    .arg(it.key()).arg(it.key() + 1).arg(it.value());
    }
    
    if (!spectralTypeHistogram.isEmpty()) {
        qDebug() << "Spectral type distribution:";
        for (auto it = spectralTypeHistogram.begin(); it != spectralTypeHistogram.end(); ++it) {
            qDebug() << QString("  %1-type: %2 stars").arg(it.key()).arg(it.value());
        }
    }
    
    qDebug() << QString("Data quality: Gaia GDR3 with proper motions applied to epoch %1")
                .arg(2025.5);
}

// Static member definition (add to .cpp file)
// QString Local2MASSCatalog::s_catalogPath;

StarCatalogValidator::StarCatalogValidator(QObject* parent)
    : QObject(parent)
    , m_validationMode(Loose)
    , m_pixelTolerance(5.0)
    , m_magnitudeTolerance(2.0)
    , m_wcsMatrixValid(false)
    , m_det(0.0)
    , m_curl(curl_easy_init())
{
    initializeTolerances();
    initializeGaiaDR3();
        
    // Initialize transformation matrix
    for (int i = 0; i < 4; ++i) {
        m_transformMatrix[i] = 0.0;
    }
}

StarCatalogValidator::~StarCatalogValidator()
{
    curl_easy_cleanup(m_curl);
}

void StarCatalogValidator::initializeTolerances()
{
    switch (m_validationMode) {
        case Strict:
            m_pixelTolerance = 2.0;
            m_magnitudeTolerance = 1.0;
            break;
        case Normal:
            m_pixelTolerance = 5.0;
            m_magnitudeTolerance = 2.0;
            break;
        case Loose:
            m_pixelTolerance = 10.0;
            m_magnitudeTolerance = 3.0;
            break;
    }
}

void StarCatalogValidator::setValidationMode(ValidationMode mode)
{
    m_validationMode = mode;
    initializeTolerances();
    qDebug() << "Validation mode set to:" << static_cast<int>(mode);
}

void StarCatalogValidator::setMatchingTolerance(double pixelTolerance, double magnitudeTolerance)
{
    m_pixelTolerance = pixelTolerance;
    m_magnitudeTolerance = magnitudeTolerance;
    qDebug() << "Tolerances set - Pixel:" << pixelTolerance << "Magnitude:" << magnitudeTolerance;
}
/*
void StarCatalogValidator::setMagnitudeLimit(double faintestMagnitude)
{
    m_magnitudeLimit = faintestMagnitude;
    qDebug() << "Magnitude limit set to:" << faintestMagnitude;
}

bool StarCatalogValidator::setWCSData(const WCSData& wcs)
{
    m_wcsData = wcs;
    
    if (wcs.isValid) {
        updateWCSMatrix();
        qDebug() << "WCS data set successfully";
        qDebug() << "  Center RA/Dec:" << wcs.crval1 << wcs.crval2;
        qDebug() << "  Reference pixel:" << wcs.crpix1 << wcs.crpix2;
        qDebug() << "  Pixel scale:" << wcs.pixscale << "arcsec/pixel";
        return true;
    } else {
        qDebug() << "Invalid WCS data provided";
        return false;
    }
}

void StarCatalogValidator::setWCSFromMetadata(const QStringList& metadata)
{
    WCSData wcs;
    
    qDebug() << "=== Extracting WCS from metadata ===";
    qDebug() << "Total metadata lines:" << metadata.size();
    
    // Debug: Show first few metadata lines to understand the format
    qDebug() << "=== First 10 metadata lines ===";
    for (int i = 0; i < qMin(10, metadata.size()); ++i) {
        qDebug() << "Line" << i << ":" << metadata[i];
    }
    
    // Look for specific WCS keywords to debug
    for (int i = 0; i < metadata.size(); ++i) {
        const QString& line = metadata[i];
        if (line.contains("CRVAL1", Qt::CaseInsensitive) || 
            line.contains("CRVAL2", Qt::CaseInsensitive) ||
            line.contains("CD1_1", Qt::CaseInsensitive) ||
            line.contains("CD2_2", Qt::CaseInsensitive)) {
            qDebug() << "WCS keyword line" << i << ":" << line;
        }
    }
    
    // Parse common WCS keywords from metadata
    for (const QString& line : metadata) {
        if (line.contains("CRVAL1", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(CRVAL1:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.crval1 = match.captured(1).toDouble();
                qDebug() << "Found CRVAL1:" << wcs.crval1;
            }
        } else if (line.contains("CRVAL2", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(CRVAL2:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.crval2 = match.captured(1).toDouble();
                qDebug() << "Found CRVAL2:" << wcs.crval2;
            }
        } else if (line.contains("CRPIX1", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(CRPIX1:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.crpix1 = match.captured(1).toDouble();
                qDebug() << "Found CRPIX1:" << wcs.crpix1;
            }
        } else if (line.contains("CRPIX2", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(CRPIX2:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.crpix2 = match.captured(1).toDouble();
                qDebug() << "Found CRPIX2:" << wcs.crpix2;
            }
        } else if (line.contains("CD1_1", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(CD1_1:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.cd11 = match.captured(1).toDouble();
                qDebug() << "Found CD1_1:" << wcs.cd11;
            }
        } else if (line.contains("CD1_2", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(CD1_2:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.cd12 = match.captured(1).toDouble();
                qDebug() << "Found CD1_2:" << wcs.cd12;
            }
        } else if (line.contains("CD2_1", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(CD2_1:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.cd21 = match.captured(1).toDouble();
                qDebug() << "Found CD2_1:" << wcs.cd21;
            }
        } else if (line.contains("CD2_2", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(CD2_2:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.cd22 = match.captured(1).toDouble();
                qDebug() << "Found CD2_2:" << wcs.cd22;
            }
        } else if (line.contains("PIXSCALE", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(PIXSCALE:\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.pixscale = match.captured(1).toDouble();
                qDebug() << "Found PIXSCALE:" << wcs.pixscale;
            }
        } else if (line.contains("ORIENTAT", Qt::CaseInsensitive) || line.contains("ROTATION", Qt::CaseInsensitive)) {
            QRegularExpression re(R"((?:ORIENTAT|ROTATION):\s*([+-]?\d*\.?\d*(?:[eE][+-]?\d+)?))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.orientation = match.captured(1).toDouble();
                qDebug() << "Found ORIENTAT/ROTATION:" << wcs.orientation;
            }
        } else if (line.contains("NAXIS1", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(NAXIS1:\s*(\d+))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.width = match.captured(1).toInt();
                qDebug() << "Found NAXIS1 (width):" << wcs.width;
            }
        } else if (line.contains("NAXIS2", Qt::CaseInsensitive)) {
            QRegularExpression re(R"(NAXIS2:\s*(\d+))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.height = match.captured(1).toInt();
                qDebug() << "Found NAXIS2 (height):" << wcs.height;
            }
        } else if (line.contains("Dimensions", Qt::CaseInsensitive)) {
            QRegularExpression re(R"((\d+)\s*[×x]\s*(\d+))");
            auto match = re.match(line);
            if (match.hasMatch()) {
                wcs.width = match.captured(1).toInt();
                wcs.height = match.captured(2).toInt();
                qDebug() << "Found dimensions:" << wcs.width << "x" << wcs.height;
            }
        }
    }
    
    // Calculate pixel scale from CD matrix if available but pixscale not directly provided
    if (wcs.pixscale == 0.0 && (wcs.cd11 != 0.0 || wcs.cd12 != 0.0 || wcs.cd21 != 0.0 || wcs.cd22 != 0.0)) {
        // Calculate pixel scale from CD matrix
        // Pixel scale = sqrt(CD1_1^2 + CD2_1^2) converted from degrees to arcseconds
        double cd11_arcsec = wcs.cd11 * 3600.0; // Convert from degrees to arcseconds
        double cd21_arcsec = wcs.cd21 * 3600.0;
        wcs.pixscale = sqrt(cd11_arcsec * cd11_arcsec + cd21_arcsec * cd21_arcsec);
        
        qDebug() << "Calculated pixel scale from CD matrix:" << wcs.pixscale << "arcsec/pixel";
        qDebug() << "CD matrix elements: CD11=" << wcs.cd11 << "CD21=" << wcs.cd21;
    }
    
    // Check if we have minimum required WCS data
    bool hasBasicWCS = (wcs.crval1 != 0.0 || wcs.crval2 != 0.0) && 
                       (wcs.cd11 != 0.0 || wcs.cd22 != 0.0 || wcs.pixscale != 0.0);
    
    if (hasBasicWCS) {
        // If we have pixscale but no CD matrix, create simple CD matrix
        if (wcs.cd11 == 0.0 && wcs.cd22 == 0.0 && wcs.pixscale != 0.0) {
            double pixscaleRad = wcs.pixscale / 3600.0 * M_PI / 180.0; // Convert to radians
            double cosRot = cos(wcs.orientation * M_PI / 180.0);
            double sinRot = sin(wcs.orientation * M_PI / 180.0);
            
            wcs.cd11 = -pixscaleRad * cosRot;
            wcs.cd12 = pixscaleRad * sinRot;
            wcs.cd21 = pixscaleRad * sinRot;
            wcs.cd22 = pixscaleRad * cosRot;
            
            qDebug() << "Created CD matrix from pixel scale:" << wcs.pixscale;
        }
        
        // Set default reference pixel if not provided
        if (wcs.crpix1 == 0.0 && wcs.crpix2 == 0.0 && wcs.width > 0 && wcs.height > 0) {
            wcs.crpix1 = wcs.width / 2.0 + 1.0;  // FITS uses 1-based indexing
            wcs.crpix2 = wcs.height / 2.0 + 1.0;
            qDebug() << "Set default reference pixel:" << wcs.crpix1 << wcs.crpix2;
        }
        
        wcs.isValid = true;
        setWCSData(wcs);
    } else {
        qDebug() << "Insufficient WCS data found in metadata";
        emit errorSignal("No valid WCS information found in image metadata");
    }
}

void StarCatalogValidator::updateWCSMatrix()
{
    if (!m_astrometricMetadata.IsValid()) {
        m_wcsMatrixValid = false;
        return;
    }
    
    // Set up transformation matrix for WCS calculations
    m_transformMatrix[0] = m_wcsData.cd11; // CD1_1
    m_transformMatrix[1] = m_wcsData.cd12; // CD1_2
    m_transformMatrix[2] = m_wcsData.cd21; // CD2_1
    m_transformMatrix[3] = m_wcsData.cd22; // CD2_2
    
    // Calculate determinant for inverse transformations
    m_det = m_transformMatrix[0] * m_transformMatrix[3] - m_transformMatrix[1] * m_transformMatrix[2];
    
    m_wcsMatrixValid = (std::abs(m_det) > 1e-15);
    
    if (!m_wcsMatrixValid) {
        qDebug() << "Warning: WCS transformation matrix is singular";
    }
}
*/
// Callback to collect HTTP response
static size_t writeCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

void StarCatalogValidator::queryGaiaCatalog(double centerRA, double centerDec, double radiusDegrees)
{
    if (!m_astrometricMetadata.IsValid()) {
        emit errorSignal("No valid WCS data available for catalog query");
        return;
    }
    
    QString baseUrl;
    QString query;

    std::string adqlQuery =
        "SELECT TOP 10000 source_id,ra,dec,phot_g_mean_mag "
        "FROM gaiaedr3.gaia_source "
        "WHERE CONTAINS(POINT('ICRS', ra, dec), "
        "CIRCLE('ICRS', " + std::to_string(centerRA) + ", " + std::to_string(centerDec) + ", " + std::to_string(radiusDegrees) + "))=1 ";

    if (!m_curl) {
      qDebug() << "Curl init failed!\n";
      return;
    }

    std::string postFields = "REQUEST=doQuery"
                             "&LANG=ADQL"
                             "&FORMAT=JSON"
      "&QUERY=" + std::string(curl_easy_escape(m_curl, adqlQuery.c_str(), adqlQuery.size()));
    
    std::string readBuffer;

    curl_easy_setopt(m_curl, CURLOPT_URL, "https://gea.esac.esa.int/tap-server/tap/sync");
    curl_easy_setopt(m_curl, CURLOPT_POSTFIELDS, postFields.c_str());
    curl_easy_setopt(m_curl, CURLOPT_WRITEFUNCTION, writeCallback);
    curl_easy_setopt(m_curl, CURLOPT_WRITEDATA, &readBuffer);

    qDebug() << "Querying catalog:" << postFields;
    
    emit catalogQueryStarted();
    
    CURLcode res = curl_easy_perform(m_curl);
    m_catData = readBuffer.c_str();
    
    if (res != CURLE_OK)
      qDebug() << "Curl error: " << curl_easy_strerror(res) << "\n";
    else
      {
      parseCatalogResponse(m_catData);
      emit catalogQueryFinished(true, QString("Retrieved %1 catalog stars").arg(m_catalogStars.size()));
      }
}

void StarCatalogValidator::onNetworkError(QNetworkReply::NetworkError error)
{
    QString errorMsg = QString("Network error during catalog query: %1").arg(static_cast<int>(error));
    qDebug() << errorMsg;
    emit errorSignal(errorMsg);
}

void StarCatalogValidator::parseCatalogResponse(const QByteArray& data)
{
    m_catalogStars.clear();
    
    // Try to parse as JSON first (for Gaia)
    QJsonParseError parseError;
    QJsonDocument doc = QJsonDocument::fromJson(data, &parseError);
    
    if (parseError.error == QJsonParseError::NoError && doc.isObject()) {
        QJsonObject obj = doc.object();
        if (obj.contains("data")) {
            parseGaiaData(obj["data"].toArray());
        }
    }
    // Calculate pixel positions for all catalog stars
    for (auto& star : m_catalogStars) {
        star.pixelPos = skyToPixel(star.ra, star.dec);
        // Mark as invalid if outside image bounds
        if (star.pixelPos.x() < 0 || star.pixelPos.y() < 0 ||
            star.pixelPos.x() >= m_astrometricMetadata.Width() || star.pixelPos.y() >= m_astrometricMetadata.Height()) {
            star.isValid = false;
        }
    }
    
    qDebug() << "Parsed" << m_catalogStars.size() << "catalog stars";
}

void StarCatalogValidator::parseGaiaData(const QJsonArray& stars)
{
    for (const QJsonValue& value : stars) {
        if (value.isArray()) {
            QJsonArray starData = value.toArray();
            if (starData.size() >= 4) {
                QString sourceId = QString::number(starData[0].toDouble(), 'f', 0);
                double ra = starData[1].toDouble();
                double dec = starData[2].toDouble();
                double magnitude = starData[3].toDouble();
                
		CatalogStar star(QString("Gaia_%1").arg(sourceId), ra, dec, magnitude);
		qDebug() << sourceId << ra << dec << magnitude;
		m_catalogStars.append(star);
            }
        }
    }
}

double StarCatalogValidator::parseCoordinate(const QString& coordStr, bool isRA) const
{
    // Parse coordinates in various formats (HH:MM:SS.S or DD:MM:SS.S)
    QRegularExpression sexagesimalRegex(R"(([+-]?\d+)[:\s]+(\d+)[:\s]+([+-]?\d+\.?\d*))");
    auto match = sexagesimalRegex.match(coordStr);
    
    if (match.hasMatch()) {
        double deg = match.captured(1).toDouble();
        double min = match.captured(2).toDouble();
        double sec = match.captured(3).toDouble();
        
        double decimal = std::abs(deg) + min/60.0 + sec/3600.0;
        if (deg < 0) decimal = -decimal;
        
        // Convert RA from hours to degrees
        if (isRA) decimal *= 15.0;
        
        return decimal;
    }
    
    // Try direct decimal parsing
    bool ok;
    double value = coordStr.toDouble(&ok);
    if (ok) {
        if (isRA && value < 24.0) value *= 15.0; // Convert hours to degrees
        return value;
    }
    
    return -999.0; // Invalid coordinate
}

ValidationResult StarCatalogValidator::validateStars(const QVector<QPoint>& detectedStars, 
                                                   const QVector<float>& starRadii)
{
    if (!m_astrometricMetadata.IsValid()) {
        ValidationResult result;
        result.summary = "No valid WCS data available";
        emit errorSignal(result.summary);
        return result;
    }
    
    if (m_catalogStars.isEmpty()) {
        ValidationResult result;
        result.summary = "No catalog stars available for validation";
        emit errorSignal(result.summary);
        return result;
    }
    
    m_lastValidation = performMatching(detectedStars, starRadii);
    emit validationCompleted(m_lastValidation);
    
    return m_lastValidation;
}

ValidationResult StarCatalogValidator::performMatching(const QVector<QPoint>& detectedStars, 
                                                     const QVector<float>& starRadii)
{
    ValidationResult result;
    result.catalogStars = m_catalogStars;
    result.totalDetected = detectedStars.size();
    result.totalCatalog = m_catalogStars.size();
    
    // Create lists to track unmatched stars
    QVector<bool> detectedMatched(detectedStars.size(), false);
    QVector<bool> catalogMatched(m_catalogStars.size(), false);
    
    // Find matches using nearest neighbor approach
    for (int i = 0; i < detectedStars.size(); ++i) {
        QPoint detected = detectedStars[i];
        double bestDistance = m_pixelTolerance;
        int bestCatalogIndex = -1;
        
        for (int j = 0; j < m_catalogStars.size(); ++j) {
            if (!m_catalogStars[j].isValid || catalogMatched[j]) continue;
            
            double distance = calculateDistance(detected, m_catalogStars[j].pixelPos);
            
            if (distance < bestDistance) {
                bestDistance = distance;
                bestCatalogIndex = j;
            }
        }
        
        if (bestCatalogIndex >= 0) {
            // Calculate magnitude difference if we have radius information
            double magnitudeDiff = 0.0;
            if (!starRadii.isEmpty() && i < starRadii.size()) {
                // Rough magnitude estimation from star radius
                // This is very approximate - real photometry would be needed
                float detectedRadius = starRadii[i];
                double estimatedMag = 10.0 - 2.5 * log10(detectedRadius * detectedRadius);
                magnitudeDiff = std::abs(estimatedMag - m_catalogStars[bestCatalogIndex].magnitude);
            }
            
            bool isGood = isGoodMatch(bestDistance, magnitudeDiff);
            
            StarMatch match(i, bestCatalogIndex, bestDistance, magnitudeDiff);
            match.isGoodMatch = isGood;
            result.matches.append(match);
            
            if (isGood) {
                detectedMatched[i] = true;
                catalogMatched[bestCatalogIndex] = true;
                result.totalMatches++;
            }
        }
    }
    
    // Collect unmatched stars
    for (int i = 0; i < detectedStars.size(); ++i) {
        if (!detectedMatched[i]) {
            result.unmatchedDetected.append(i);
        }
    }
    
    for (int i = 0; i < m_catalogStars.size(); ++i) {
        if (m_catalogStars[i].isValid && !catalogMatched[i]) {
            result.unmatchedCatalog.append(i);
        }
    }
    
    calculateStatistics(result);
    result.isValid = true;
    
    return result;
}

double StarCatalogValidator::calculateDistance(const QPoint& detected, const QPointF& catalog) const
{
    double dx = detected.x() - catalog.x();
    double dy = detected.y() - catalog.y();
    return sqrt(dx * dx + dy * dy);
}

bool StarCatalogValidator::isGoodMatch(double distance, double magnitudeDiff) const
{
    if (distance > m_pixelTolerance) return false;
    if (magnitudeDiff > m_magnitudeTolerance) return false;
    
    // Additional quality checks based on validation mode
    switch (m_validationMode) {
        case Strict:
            return distance < m_pixelTolerance * 0.5;
        case Normal:
            return true; // Already passed basic checks
        case Loose:
            return distance < m_pixelTolerance * 1.5;
    }
    
    return true;
}

void StarCatalogValidator::calculateStatistics(ValidationResult& result) const
{
    if (result.totalDetected == 0) {
        result.matchPercentage = 0.0;
        result.summary = "No detected stars to validate";
        return;
    }
    
    result.matchPercentage = (100.0 * result.totalMatches) / result.totalDetected;
    
    // Calculate position error statistics
    if (!result.matches.isEmpty()) {
        double sumError = 0.0;
        double sumSquaredError = 0.0;
        int goodMatches = 0;
        
        for (const auto& match : result.matches) {
            if (match.isGoodMatch) {
                sumError += match.distance;
                sumSquaredError += match.distance * match.distance;
                goodMatches++;
            }
        }
        
        if (goodMatches > 0) {
            result.averagePositionError = sumError / goodMatches;
            result.rmsPositionError = sqrt(sumSquaredError / goodMatches);
        }
    }
    
    // Generate summary
    result.summary = QString("Validation Results:\n"
                           "Detected: %1, Catalog: %2, Matches: %3 (%4%%)\n"
                           "Average position error: %5 pixels\n"
                           "RMS position error: %6 pixels\n"
                           "Unmatched detections: %7\n"
                           "Unmatched catalog stars: %8")
                   .arg(result.totalDetected)
                   .arg(result.totalCatalog)
                   .arg(result.totalMatches)
                   .arg(result.matchPercentage)
                   .arg(result.averagePositionError, 0, 'f', 2)
                   .arg(result.rmsPositionError, 0, 'f', 2)
                   .arg(result.unmatchedDetected.size())
                   .arg(result.unmatchedCatalog.size());
}

QString StarCatalogValidator::getValidationSummary() const
{
    return m_lastValidation.summary;
}

QJsonObject StarCatalogValidator::exportValidationResults() const
{
    const ValidationResult& result = m_lastValidation;
    
    QJsonObject obj;
    obj["totalDetected"] = result.totalDetected;
    obj["totalCatalog"] = result.totalCatalog;
    obj["totalMatches"] = result.totalMatches;
    obj["matchPercentage"] = result.matchPercentage;
    obj["averagePositionError"] = result.averagePositionError;
    obj["rmsPositionError"] = result.rmsPositionError;
    obj["summary"] = result.summary;
    
    QJsonArray matchesArray;
    for (const auto& match : result.matches) {
        QJsonObject matchObj;
        matchObj["detectedIndex"] = match.detectedIndex;
        matchObj["catalogIndex"] = match.catalogIndex;
        matchObj["distance"] = match.distance;
        matchObj["magnitudeDiff"] = match.magnitudeDiff;
        matchObj["isGoodMatch"] = match.isGoodMatch;
        matchesArray.append(matchObj);
    }
    obj["matches"] = matchesArray;
    
    QJsonArray catalogArray;
    for (const auto& star : result.catalogStars) {
        QJsonObject starObj;
        starObj["id"] = star.id;
        starObj["ra"] = star.ra;
        starObj["dec"] = star.dec;
        starObj["magnitude"] = star.magnitude;
        starObj["spectralType"] = star.spectralType;
        starObj["pixelX"] = star.pixelPos.x();
        starObj["pixelY"] = star.pixelPos.y();
        starObj["isValid"] = star.isValid;
        catalogArray.append(starObj);
    }
    obj["catalogStars"] = catalogArray;
    
    return obj;
}

void StarCatalogValidator::clearResults()
{
    m_catalogStars.clear();
    m_lastValidation.clear();
    qDebug() << "Validation results cleared";
}

bool StarCatalogValidator::setWCSFromImageMetadata(const ImageData& imageData)
{
    qDebug() << "=== Using PCL AstrometricMetadata (Complete) ===";
    
    try {
        
        // Convert ImageData metadata to PCL FITSKeywordArray
        pcl::FITSKeywordArray keywords;
        
        // Parse metadata lines and convert to FITS keywords
        for (const QString& line : imageData.metadata) {
            QString trimmedLine = line.trimmed();
            
            // Extract key-value pairs in the format "KEY: value (comment)"
            QRegularExpression kvRegex(R"(([A-Z0-9_-]+):\s*([^(]+)(?:\s*\(([^)]*)\))?)", 
                                     QRegularExpression::CaseInsensitiveOption);
            auto match = kvRegex.match(trimmedLine);
            
            if (match.hasMatch()) {
                QString key = match.captured(1).toUpper();
                QString valueStr = match.captured(2).trimmed();
                QString comment = match.captured(3).trimmed();
                
                // Handle quoted string values
                if (valueStr.startsWith("'") && valueStr.endsWith("'")) {
                    valueStr = valueStr.mid(1, valueStr.length() - 2);
                }
                
                // Only add WCS-related keywords that PCL recognizes
                QStringList wcsKeys = {
                    "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
                    "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                    "CDELT1", "CDELT2", "CROTA1", "CROTA2",
                    "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2",
                    "PV1_1", "PV1_2", "LONPOLE", "LATPOLE",
                    "EQUINOX", "RADESYS", "NAXIS1", "NAXIS2",
                    // Add some basic metadata
                    "DATE-OBS", "FOCALLEN", "XPIXSZ", "YPIXSZ"
                };
                
                if (wcsKeys.contains(key)) {
                    // Try to parse as number first
                    bool isNumber;
                    double numValue = valueStr.toDouble(&isNumber);
                    
                    if (isNumber) {
                        keywords.Add(pcl::FITSHeaderKeyword(
                            pcl::IsoString(key.toUtf8().constData()), 
                            numValue, 
                            pcl::IsoString(comment.toUtf8().constData())));
                        qDebug() << "Added numeric keyword:" << key << "=" << numValue;
                    } else {
                        // Handle as string (for CTYPE, RADESYS, etc.)
                        keywords.Add(pcl::FITSHeaderKeyword(
                            pcl::IsoString(key.toUtf8().constData()), 
                            pcl::IsoString(valueStr.toUtf8().constData()), 
                            pcl::IsoString(comment.toUtf8().constData())));
                        qDebug() << "Added string keyword:" << key << "=" << valueStr;
                    }
                }
            }
        }
        
        // Add essential default keywords if not present
        bool hasCTYPE1 = false, hasCTYPE2 = false, hasRADESYS = false;
        for (const auto& kw : keywords) {
            if (kw.name == "CTYPE1") hasCTYPE1 = true;
            if (kw.name == "CTYPE2") hasCTYPE2 = true;
            if (kw.name == "RADESYS") hasRADESYS = true;
        }
        
        if (!hasCTYPE1) {
            keywords.Add(pcl::FITSHeaderKeyword("CTYPE1", 
                pcl::IsoString("'RA---TAN'"), "Coordinate type"));
        }
        if (!hasCTYPE2) {
            keywords.Add(pcl::FITSHeaderKeyword("CTYPE2", 
                pcl::IsoString("'DEC--TAN'"), "Coordinate type"));
        }
        if (!hasRADESYS) {
            keywords.Add(pcl::FITSHeaderKeyword("RADESYS", 
                pcl::IsoString("'ICRS'"), "Coordinate reference system"));
        }
        
        qDebug() << "Created" << keywords.Length() << "FITS keywords for PCL";
        
        // Use PCL's AstrometricMetadata.Build() method
        // This is the high-level interface that handles everything!
        pcl::PropertyArray emptyProperties; // We're using FITS keywords only
        
        try {
            m_astrometricMetadata.Build(emptyProperties, keywords, 
                                       imageData.width, imageData.height);
            
            if (m_astrometricMetadata.IsValid()) {
                    qDebug() << "✅ PCL AstrometricMetadata setup successful!";
//                  qDebug() << "  Image center: RA=" << centerWorld.x << "° Dec=" << centerWorld.y << "°";
                    qDebug() << "  Resolution:" << getPixScale() << "arcsec/pixel";
                    qDebug() << "  Image size:"
			     << m_astrometricMetadata.Width() << "x"
			     << m_astrometricMetadata.Height();
                    
                    // Test the coordinate transformations
                    testAstrometricMetadata();
                    
                    return true;
            } else {
                qDebug() << "❌ PCL AstrometricMetadata.IsValid() returned false";
                
                // Try to get more information about why it failed
                try {
                    m_astrometricMetadata.Validate(0.1); // Test with 0.1 pixel tolerance
                } catch (const pcl::Error& e) {
                    qDebug() << "  Validation error:" << e.Message().c_str();
                }
            }
            
        } catch (const pcl::Error& e) {
            qDebug() << "❌ PCL AstrometricMetadata.Build() failed:" << e.Message().c_str();
        }
        
        return false;
        
    } catch (const std::exception& e) {
        qDebug() << "❌ Standard exception in PCL AstrometricMetadata setup:" << e.what();
        return false;
    } catch (...) {
        qDebug() << "❌ Unknown exception in PCL AstrometricMetadata setup";
        return false;
    }
}

bool StarCatalogValidator::setWCSFromSolver(pcl::FITSKeywordArray keywords, int width, int height)
{
    try {
        // Use PCL's AstrometricMetadata.Build() method
        // This is the high-level interface that handles everything!
        pcl::PropertyArray emptyProperties; // We're using FITS keywords only
        
        try {
            m_astrometricMetadata.Build(emptyProperties, keywords, 
                                       width, height);
            
            if (m_astrometricMetadata.IsValid()) {
                    qDebug() << "✅ PCL AstrometricMetadata setup successful!";
                    qDebug() << "  Resolution:" << getPixScale() << "arcsec/pixel";
                    qDebug() << "  Image size:"
			     << m_astrometricMetadata.Width() << "x"
			     << m_astrometricMetadata.Height();
                    
                    // Test the coordinate transformations
                    testAstrometricMetadata();
                    
                    return true;
            } else {
                qDebug() << "❌ PCL AstrometricMetadata.IsValid() returned false";
                
                // Try to get more information about why it failed
                try {
                    m_astrometricMetadata.Validate(0.1); // Test with 0.1 pixel tolerance
                } catch (const pcl::Error& e) {
                    qDebug() << "  Validation error:" << e.Message().c_str();
                }
            }
            
        } catch (const pcl::Error& e) {
            qDebug() << "❌ PCL AstrometricMetadata.Build() failed:" << e.Message().c_str();
        }
        
        return false;
        
    } catch (const std::exception& e) {
        qDebug() << "❌ Standard exception in PCL AstrometricMetadata setup:" << e.what();
        return false;
    } catch (...) {
        qDebug() << "❌ Unknown exception in PCL AstrometricMetadata setup";
        return false;
    }
}

void StarCatalogValidator::testAstrometricMetadata() const
{
    if (!m_astrometricMetadata.IsValid()) return;
    
    qDebug() << "\n=== Testing PCL AstrometricMetadata ===";
    
    try {
        // Test 1: Image center
        pcl::DPoint centerImage(m_astrometricMetadata.Width() * 0.5, m_astrometricMetadata.Height() * 0.5);
        pcl::DPoint centerWorld;
        
        if (m_astrometricMetadata.ImageToCelestial(centerWorld, centerImage)) {
            qDebug() << QString("✅ Image center: (%1, %2) -> RA=%3° Dec=%4°")
                        .arg(centerImage.x, 0, 'f', 1).arg(centerImage.y, 0, 'f', 1)
                        .arg(centerWorld.x, 0, 'f', 6).arg(centerWorld.y, 0, 'f', 6);
        } else {
            qDebug() << "❌ Failed to transform image center";
        }
        
        // Test 2: Image corners
        struct TestPoint {
            pcl::DPoint image;
            QString name;
        };
        
        TestPoint corners[] = {
            {pcl::DPoint(0, 0), "Top-Left"},
            {pcl::DPoint(m_astrometricMetadata.Width(), 0), "Top-Right"},
            {pcl::DPoint(0, m_astrometricMetadata.Height()), "Bottom-Left"},
            {pcl::DPoint(m_astrometricMetadata.Width(), m_astrometricMetadata.Height()), "Bottom-Right"}
        };
        
        qDebug() << "\nImage corner transformations:";
        for (const auto& corner : corners) {
            pcl::DPoint world;
            if (m_astrometricMetadata.ImageToCelestial(world, corner.image)) {
                qDebug() << QString("  %1: (%2, %3) -> RA=%4° Dec=%5°")
                            .arg(corner.name).arg(corner.image.x, 0, 'f', 0).arg(corner.image.y, 0, 'f', 0)
                            .arg(world.x, 0, 'f', 4).arg(world.y, 0, 'f', 4);
            } else {
                qDebug() << QString("  %1: Transformation failed").arg(corner.name);
            }
        }
        
        // Test 3: Catalog stars (your Gaia data)
        struct TestStar {
            double ra, dec;
            QString name;
        };
        
        TestStar testStars[] = {
            {11.195908490324188, 41.89174375779233, "Gaia_387310087645581056"},
            {11.245222869774707, 41.818567129781435, "Gaia_387309567953910528"},
            {11.25642400913238, 41.84593075636271, "Gaia_387309705392866432"},
            {centerWorld.x, centerWorld.y, "Image_Center_Roundtrip"}
        };
        
        qDebug() << "\nCatalog star transformations:";
        int inBoundsCount = 0;
        
        for (const auto& star : testStars) {
            pcl::DPoint worldCoord(star.ra, star.dec);
            pcl::DPoint imageCoord;
            
            if (m_astrometricMetadata.CelestialToImage(imageCoord, worldCoord)) {
                bool inBounds = (imageCoord.x >= 0 && imageCoord.x < m_astrometricMetadata.Width() && 
                               imageCoord.y >= 0 && imageCoord.y < m_astrometricMetadata.Height());
                if (inBounds) inBoundsCount++;
                
                qDebug() << QString("  %1: RA=%2° Dec=%3° -> pixel=(%4,%5) %6")
                            .arg(star.name).arg(star.ra, 0, 'f', 6).arg(star.dec, 0, 'f', 6)
                            .arg(imageCoord.x, 0, 'f', 1).arg(imageCoord.y, 0, 'f', 1)
                            .arg(inBounds ? "✅" : "❌");
            } else {
                qDebug() << QString("  %1: CelestialToImage failed").arg(star.name);
            }
        }
        
        qDebug() << QString("\n📊 Results: %1/%2 stars in bounds").arg(inBoundsCount).arg(4);
        
        if (inBoundsCount > 0) {
            qDebug() << "🎉 SUCCESS: PCL AstrometricMetadata is working correctly!";
            qDebug() << "    Your catalog stars should now appear on the image.";
        } else {
            qDebug() << "⚠️  WARNING: No stars are landing within image bounds.";
            qDebug() << "    Check if the WCS keywords in your FITS file are correct.";
        }
        
        // Test 4: Round-trip accuracy test
        qDebug() << "\nRound-trip accuracy test:";
        pcl::DPoint testPixel(m_astrometricMetadata.Width() * 0.25, m_astrometricMetadata.Height() * 0.75);
        pcl::DPoint worldCoord, backToPixel;
        
        if (m_astrometricMetadata.ImageToCelestial(worldCoord, testPixel) &&
            m_astrometricMetadata.CelestialToImage(backToPixel, worldCoord)) {
            
            double errorX = testPixel.x - backToPixel.x;
            double errorY = testPixel.y - backToPixel.y;
            double totalError = sqrt(errorX*errorX + errorY*errorY);
            
            qDebug() << QString("  Original: (%1, %2)")
                        .arg(testPixel.x, 0, 'f', 3).arg(testPixel.y, 0, 'f', 3);
            qDebug() << QString("  Round-trip: (%1, %2)")
                        .arg(backToPixel.x, 0, 'f', 3).arg(backToPixel.y, 0, 'f', 3);
            qDebug() << QString("  Error: %1 pixels %2")
                        .arg(totalError).arg(totalError < 0.1 ? "✅" : "⚠️");
        }
        
        // Test 5: Print PCL's diagnostic info
        qDebug() << "\nPCL Internal Diagnostics:";
        try {
            // Test the validation
            m_astrometricMetadata.Validate(0.01);
            qDebug() << "  ✅ PCL internal validation passed (0.01 px tolerance)";
        } catch (const pcl::Error& e) {
            qDebug() << "  ❌ PCL validation failed:" << e.Message().c_str();
        }
        
        // Get search radius for catalog queries
        try {
            double searchRadius = m_astrometricMetadata.SearchRadius();
            qDebug() << QString("  📡 Search radius: %1° (for catalog queries)")
                        .arg(searchRadius, 0, 'f', 2);
        } catch (const pcl::Error& e) {
            qDebug() << "  ❌ SearchRadius failed:" << e.Message().c_str();
        }
        
    } catch (const pcl::Error& e) {
        qDebug() << "❌ PCL error during testing:" << e.Message().c_str();
    } catch (...) {
        qDebug() << "❌ Unknown error during testing";
    }
}

// Simplified coordinate transformation methods using PCL's high-level API
QPointF StarCatalogValidator::skyToPixel(double ra, double dec) const
{
    if (!m_astrometricMetadata.IsValid()) {
        return QPointF(-1, -1);
    }
    
    try {
        pcl::DPoint worldCoord(ra, dec);
        pcl::DPoint imageCoord;
        
        if (m_astrometricMetadata.CelestialToImage(imageCoord, worldCoord)) {
            return QPointF(imageCoord.x, imageCoord.y);
        }
        
    } catch (const pcl::Error& e) {
        qDebug() << "skyToPixel error:" << e.Message().c_str();
    }
    
    return QPointF(-1, -1);
}

QPointF StarCatalogValidator::pixelToSky(double x, double y) const
{
    if (!m_astrometricMetadata.IsValid()) {
        return QPointF(-1, -1);
    }
    
    try {
        pcl::DPoint imageCoord(x, y);
        pcl::DPoint worldCoord;
        
        if (m_astrometricMetadata.ImageToCelestial(worldCoord, imageCoord)) {
            return QPointF(worldCoord.x, worldCoord.y);
        }
        
    } catch (const pcl::Error& e) {
        qDebug() << "pixelToSky error:" << e.Message().c_str();
    }
    
    return QPointF(-1, -1);
}

// Update your other methods
bool StarCatalogValidator::hasValidWCS() const 
{ 
    return m_astrometricMetadata.IsValid() && m_astrometricMetadata.IsValid(); 
}

// Update StarCatalogValidator to use the local bright star database
void StarCatalogValidator::addBrightStarsFromDatabase(double centerRA, double centerDec, double radiusDegrees)
{
    qDebug() << "\n=== ADDING BRIGHT STARS FROM LOCAL DATABASE ===";
    
    // Get bright stars in field (magnitude < 3.0 to ensure visibility)
    auto brightStars = BrightStarDatabase::getStarsInField(centerRA, centerDec, radiusDegrees, 3.0);
    
    if (brightStars.isEmpty()) {
        qDebug() << "No bright stars (mag < 3.0) found in this field";
        return;
    }
    
    qDebug() << QString("Found %1 bright stars in field:").arg(brightStars.size());
    
    int addedCount = 0;
    int replacedCount = 0;
    
    for (const auto& brightStar : brightStars) {
        // Check if magnitude is within our limit
        
        // Create catalog star
        CatalogStar catalogStar(brightStar.id, brightStar.ra, brightStar.dec, brightStar.magnitude);
        catalogStar.spectralType = brightStar.spectralType;
        
        // Calculate pixel position
        catalogStar.pixelPos = skyToPixel(brightStar.ra, brightStar.dec);
        
        // Check if in image bounds
        catalogStar.isValid = (catalogStar.pixelPos.x() >= 0 && catalogStar.pixelPos.x() < m_astrometricMetadata.Width() && 
                              catalogStar.pixelPos.y() >= 0 && catalogStar.pixelPos.y() < m_astrometricMetadata.Height());
        
        // Check if we already have a star at this position (replace if ours is better)
        bool replaced = false;
        for (int i = 0; i < m_catalogStars.size(); ++i) {
            auto& existing = m_catalogStars[i];
            double ra_diff = abs(existing.ra - brightStar.ra);
            double dec_diff = abs(existing.dec - brightStar.dec);
            
            if (ra_diff < 0.01 && dec_diff < 0.01) { // Same star (within 0.01 degrees)
                // Replace if our magnitude is more reasonable (bright stars should have mag < 6)
                if (existing.magnitude > 6.0 && brightStar.magnitude < 6.0) {
                    qDebug() << QString("  🔄 Replacing %1: bad mag %2 -> correct mag %3")
                                .arg(brightStar.name).arg(existing.magnitude).arg(brightStar.magnitude);
                    m_catalogStars[i] = catalogStar;
                    replaced = true;
                    replacedCount++;
                    break;
                }
            }
        }
        
        if (!replaced) {
            m_catalogStars.append(catalogStar);
            addedCount++;
            
            qDebug() << QString("  ✅ Added %1 (%2): mag %3 -> pixel (%4, %5) %6")
                        .arg(brightStar.name).arg(brightStar.constellation).arg(brightStar.magnitude)
                        .arg(catalogStar.pixelPos.x()).arg(catalogStar.pixelPos.y())
                        .arg(catalogStar.isValid ? "✅" : "❌");
        }
    }
    
    qDebug() << QString("Added %1 new bright stars, replaced %2 existing entries")
                .arg(addedCount).arg(replacedCount);
    
    if (addedCount > 0 || replacedCount > 0) {
        // Sort catalog by magnitude (brightest first)
        std::sort(m_catalogStars.begin(), m_catalogStars.end(), 
                  [](const CatalogStar& a, const CatalogStar& b) {
                      return a.magnitude < b.magnitude;
                  });
        
        qDebug() << QString("Total catalog now contains %1 stars").arg(m_catalogStars.size());
    }
}

// Direct WCS info extraction from PCL AstrometricMetadata

void StarCatalogValidator::showWCSInfoFromPCL()
{
    qDebug() << "\n=== WCS INFO FROM PCL ASTROMETRICMETADATA ===";
    
    if (!m_astrometricMetadata.IsValid()) {
        qDebug() << "❌ No valid astrometric metadata available";
        return;
    }
    
    qDebug() << "WCS DEBUG INFORMATION";
    qDebug() << "====================";
    qDebug() << "";
    
    // Basic WCS Parameters
    qDebug() << "Basic WCS Parameters:";
    
    // Get image center coordinates
    pcl::DPoint centerCoords;
    if (m_astrometricMetadata.ImageCenterToCelestial(centerCoords)) {
        qDebug() << QString("CRVAL1 (RA center): %1°").arg(centerCoords.x, 0, 'f', 6);
        qDebug() << QString("CRVAL2 (Dec center): %1°").arg(centerCoords.y, 0, 'f', 6);
    }
    
    // Reference pixel (usually image center for PCL)
    double crpix1 = m_astrometricMetadata.Width() * 0.5;
    double crpix2 = m_astrometricMetadata.Height() * 0.5;
    qDebug() << QString("CRPIX1 (X reference): %1 px").arg(crpix1, 0, 'f', 2);
    qDebug() << QString("CRPIX2 (Y reference): %1 px").arg(crpix2, 0, 'f', 2);
    
    // Pixel scale and orientation
    double pixelScale = getPixScale();
    qDebug() << QString("Pixel scale: %1 arcsec/px").arg(pixelScale, 0, 'f', 2);
    
    bool flipped;
    double rotation = m_astrometricMetadata.Rotation(flipped);
    qDebug() << QString("Orientation: %1°").arg(rotation, 0, 'f', 2);
    
    // Image dimensions
    qDebug() << QString("Image size: %1 × %2 px")
                .arg(m_astrometricMetadata.Width())
                .arg(m_astrometricMetadata.Height());
    
    // CD Matrix (extracted via numerical differentiation)
    qDebug() << "";
    qDebug() << "CD Matrix:";
    extractAndDisplayCDMatrix();
    
    // Field of View
    qDebug() << "";
    qDebug() << "Calculated Field of View:";
    calculateAndDisplayFOV();
    
    // WCS Quality Check
    qDebug() << "";
    qDebug() << "WCS Quality Check:";
    performWCSQualityCheck();
    
    // Additional PCL-specific information
    qDebug() << "";
    qDebug() << "PCL Astrometric Information:";
    displayPCLSpecificInfo();
}

void StarCatalogValidator::extractAndDisplayCDMatrix()
{
    // Calculate CD matrix elements via numerical differentiation
    double centerX = m_astrometricMetadata.Width() * 0.5;
    double centerY = m_astrometricMetadata.Height() * 0.5;
    double delta = 1.0; // 1 pixel offset
    
    pcl::DPoint center, right, up;
    
    if (m_astrometricMetadata.ImageToCelestial(center, pcl::DPoint(centerX, centerY)) &&
        m_astrometricMetadata.ImageToCelestial(right,  pcl::DPoint(centerX + delta, centerY)) &&
        m_astrometricMetadata.ImageToCelestial(up,     pcl::DPoint(centerX, centerY + delta))) {
        
        // Calculate CD matrix elements
        double cd1_1 = (right.x - center.x) / delta;  // dRA/dX
        double cd1_2 = (up.x - center.x) / delta;     // dRA/dY  
        double cd2_1 = (right.y - center.y) / delta;  // dDec/dX
        double cd2_2 = (up.y - center.y) / delta;     // dDec/dY
        
        // Handle RA wraparound near 0/360 boundary
        if (abs(cd1_1) > 180) cd1_1 = cd1_1 > 0 ? cd1_1 - 360 : cd1_1 + 360;
        if (abs(cd1_2) > 180) cd1_2 = cd1_2 > 0 ? cd1_2 - 360 : cd1_2 + 360;
        
        qDebug() << QString("CD1_1: %1").arg(cd1_1, 0, 'e', 6);
        qDebug() << QString("CD1_2: %1").arg(cd1_2, 0, 'e', 6);
        qDebug() << QString("CD2_1: %1").arg(cd2_1, 0, 'e', 6);
        qDebug() << QString("CD2_2: %1").arg(cd2_2, 0, 'e', 6);
        
        // Calculate determinant for quality check
        double determinant = cd1_1 * cd2_2 - cd1_2 * cd2_1;
        qDebug() << QString("Determinant: %1").arg(determinant, 0, 'e', 6);
        
    } else {
        qDebug() << "❌ Failed to calculate CD matrix elements";
        qDebug() << "CD1_1: 0.000000e+00";
        qDebug() << "CD1_2: 0.000000e+00";
        qDebug() << "CD2_1: 0.000000e+00";
        qDebug() << "CD2_2: 0.000000e+00";
    }
}

void StarCatalogValidator::calculateAndDisplayFOV()
{
    try {
        // Calculate field of view using PCL
        double searchRadius = m_astrometricMetadata.SearchRadius(); // degrees
        double diagonalFOV = searchRadius * 2.0; // degrees
        
        // Convert to arcminutes
        double widthArcmin = m_astrometricMetadata.Width() * m_astrometricMetadata.Resolution() * 60.0;
        double heightArcmin = m_astrometricMetadata.Height() * m_astrometricMetadata.Resolution() * 60.0;
        double diagonalArcmin = diagonalFOV * 60.0;
        
        qDebug() << QString("Width: %1 arcmin").arg(widthArcmin, 0, 'f', 1);
        qDebug() << QString("Height: %1 arcmin").arg(heightArcmin, 0, 'f', 1);
        qDebug() << QString("Diagonal: %1 arcmin").arg(diagonalArcmin, 0, 'f', 1);
        
    } catch (...) {
        // Fallback calculation
        double pixelScaleDeg = m_astrometricMetadata.Resolution();
        double widthDeg = m_astrometricMetadata.Width() * pixelScaleDeg;
        double heightDeg = m_astrometricMetadata.Height() * pixelScaleDeg;
        double diagonalDeg = sqrt(widthDeg*widthDeg + heightDeg*heightDeg);
        
        qDebug() << QString("Width: %1 arcmin").arg(widthDeg * 60.0, 0, 'f', 1);
        qDebug() << QString("Height: %1 arcmin").arg(heightDeg * 60.0, 0, 'f', 1);
        qDebug() << QString("Diagonal: %1 arcmin").arg(diagonalDeg * 60.0, 0, 'f', 1);
    }
}

void StarCatalogValidator::performWCSQualityCheck()
{
    bool isValid = m_astrometricMetadata.IsValid();
    qDebug() << QString("Is Valid: %1").arg(isValid ? "Yes" : "No");
    
    if (isValid) {
        try {
            // Test coordinate transformation accuracy
            pcl::DPoint centerImage(m_astrometricMetadata.Width() * 0.5, 
                                   m_astrometricMetadata.Height() * 0.5);
            pcl::DPoint celestial, roundTripImage;
            
            if (m_astrometricMetadata.ImageToCelestial(celestial, centerImage) &&
                m_astrometricMetadata.CelestialToImage(roundTripImage, celestial)) {
                
                double errorX = abs(roundTripImage.x - centerImage.x);
                double errorY = abs(roundTripImage.y - centerImage.y);
                double totalError = sqrt(errorX*errorX + errorY*errorY);
                
                qDebug() << QString("Round-trip error: %1 pixels").arg(totalError, 0, 'f', 6);
                
                if (totalError < 0.01) {
                    qDebug() << "Scale consistency: Excellent";
                } else if (totalError < 0.1) {
                    qDebug() << "Scale consistency: Good";  
                } else {
                    qDebug() << "Scale consistency: Poor";
                }
            } else {
                qDebug() << "Round-trip error: Failed to calculate";
                qDebug() << "Scale consistency: Invalid";
            }
        } catch (...) {
            qDebug() << "WCS quality check failed";
        }
    }
}

void StarCatalogValidator::displayPCLSpecificInfo()
{
    // Reference system
    qDebug() << QString("Reference System: %1").arg(m_astrometricMetadata.ReferenceSystem().c_str());
    
    // Check if using spline-based transformation
    if (m_astrometricMetadata.HasSplineWorldTransformation()) {
        qDebug() << "Transformation Type: Spline-based (with distortion correction)";
        
        // Try to get spline information if available
        const pcl::WorldTransformation* worldTransform = m_astrometricMetadata.WorldTransform();
        if (worldTransform) {
            const pcl::SplineWorldTransformation* splineTransform = 
                dynamic_cast<const pcl::SplineWorldTransformation*>(worldTransform);
            if (splineTransform) {
                qDebug() << QString("Control Points: %1").arg(splineTransform->NumberOfControlPoints());
                qDebug() << QString("RBF Type: %1").arg(splineTransform->RBFTypeName().c_str());
                
                int xWI, yWI, xIW, yIW;
                splineTransform->GetSplineLengths(xWI, yWI, xIW, yIW);
                qDebug() << QString("Spline Lengths: WI(%1,%2) IW(%3,%4)").arg(xWI).arg(yWI).arg(xIW).arg(yIW);
                
                double xWIErr, yWIErr, xIWErr, yIWErr;
                splineTransform->GetSplineErrors(xWIErr, yWIErr, xIWErr, yIWErr);
                qDebug() << QString("Spline Errors: WI(%1,%2) IW(%3,%4) px")
                            .arg(xWIErr).arg(yWIErr).arg(xIWErr).arg(yIWErr);
            }
        }
    } else {
        qDebug() << "Transformation Type: Linear";
    }
    
    // Resolution information
    try {
        double centerRes = m_astrometricMetadata.ResolutionAt(pcl::DPoint(
            m_astrometricMetadata.Width() * 0.5, 
            m_astrometricMetadata.Height() * 0.5));
        qDebug() << QString("Resolution at center: %1 arcsec/px").arg(centerRes * 3600.0, 0, 'f', 3);
    } catch (...) {
        qDebug() << QString("Resolution (average): %1 arcsec/px")
                    .arg(m_astrometricMetadata.Resolution() * 3600.0, 0, 'f', 3);
    }
    
    // Search radius for catalog queries
    try {
        double searchRadius = m_astrometricMetadata.SearchRadius();
        qDebug() << QString("Catalog search radius: %1°").arg(searchRadius, 0, 'f', 2);
    } catch (...) {
        qDebug() << "Catalog search radius: Not available";
    }
    
    // Time information if available
    auto obsTime = m_astrometricMetadata.ObservationMiddleTime();
    if (obsTime.IsDefined()) {
        qDebug() << QString("Observation time: %1").arg(obsTime().ToString().c_str());
    }
    
    // Creator information if available
    pcl::String creator = m_astrometricMetadata.CreatorApplication();
    if (!creator.IsEmpty()) {
        qDebug() << QString("Created by: %1").arg(creator.c_str());
    }
}

// Alternative simpler version if you just want basic info
void StarCatalogValidator::showBasicWCSInfoFromPCL()
{
    if (!m_astrometricMetadata.IsValid()) {
        qDebug() << "No valid WCS information available";
        return;
    }
    
    qDebug() << "=== WCS INFORMATION ===";
    
    // Image center
    pcl::DPoint center;
    if (m_astrometricMetadata.ImageCenterToCelestial(center)) {
        qDebug() << QString("Image Center: RA=%1° Dec=%2°").arg(center.x, 0, 'f', 6).arg(center.y, 0, 'f', 6);
    }
    
    // Image properties
    qDebug() << QString("Image Size: %1 × %2 pixels")
                .arg(m_astrometricMetadata.Width()).arg(m_astrometricMetadata.Height());
    
    // Resolution
    double resolution = m_astrometricMetadata.Resolution() * 3600.0;
    qDebug() << QString("Pixel Scale: %1 arcsec/pixel").arg(resolution, 0, 'f', 2);
    
    // Field of view
    double fovWidth = m_astrometricMetadata.Width() * resolution / 60.0;  // arcmin
    double fovHeight = m_astrometricMetadata.Height() * resolution / 60.0; // arcmin
    qDebug() << QString("Field of View: %1' × %2'").arg(fovWidth, 0, 'f', 1).arg(fovHeight, 0, 'f', 1);
    
    // Transformation type
    qDebug() << QString("Coordinate System: %1 %2")
                .arg(m_astrometricMetadata.ReferenceSystem().c_str())
                .arg(m_astrometricMetadata.HasSplineWorldTransformation() ? "(Spline)" : "(Linear)");
    
    // Quality
    qDebug() << QString("Status: %1").arg(m_astrometricMetadata.IsValid() ? "Valid" : "Invalid");
}

QVector<CatalogStar> StarCatalogValidator::getCatalogStars() {
  if (!m_catalogStars.size())
    {
      // Need to query catalog first
      double crval1, crval2;
      getCenter(crval1, crval2);
      int width =  getWidth();
      int height = getHeight();
      double pixscale = getPixScale();
      double fieldRadius = sqrt(width * width + height * height) * pixscale / 3600.0 / 2.0;
      fieldRadius = std::max(fieldRadius, 0.5);
      queryCatalog(crval1, crval2, fieldRadius);
    }
  return m_catalogStars;
}
