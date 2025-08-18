#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

const double DEG2RAD = M_PI / 180.0;
const double RAD2DEG = 180.0 / M_PI;
const double AU = 149597870.7; // km

struct CometData {
    double jd;
    double eccentricity;
    double periapsis_q;    // AU
    double inclination;    // degrees
    double long_node;      // degrees
    double arg_peri;       // degrees
    double true_anomaly;   // degrees
};

struct EarthData {
    double jd;
    double eccentricity;
    double periapsis_q;    // AU
    double inclination;    // degrees
    double long_node;      // degrees
    double arg_peri;       // degrees
    double true_anomaly;   // degrees
};

class GeocentricCalculator {
private:
    std::vector<CometData> comet_elements;
    std::vector<EarthData> earth_elements;
    
public:
    bool loadCometFile(const std::string& filename) {
        std::cout << "Loading comet file: " << filename << std::endl;
        return loadHorizonsFile(filename, comet_elements, "comet");
    }
    
    bool loadEarthFile(const std::string& filename) {
        std::cout << "Loading Earth file: " << filename << std::endl;
        std::vector<CometData> temp_earth;
        bool result = loadHorizonsFile(filename, temp_earth, "Earth");
        
        // Convert to EarthData format
        for (const auto& data : temp_earth) {
            EarthData earth;
            earth.jd = data.jd;
            earth.eccentricity = data.eccentricity;
            earth.periapsis_q = data.periapsis_q;
            earth.inclination = data.inclination;
            earth.long_node = data.long_node;
            earth.arg_peri = data.arg_peri;
            earth.true_anomaly = data.true_anomaly;
            earth_elements.push_back(earth);
        }
        
        return result;
    }
    
private:
    template<typename T>
    bool loadHorizonsFile(const std::string& filename, std::vector<T>& elements, const std::string& type) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "Cannot open file: " << filename << std::endl;
            return false;
        }
        
        std::string line;
        bool in_data = false;
        int count = 0;
        
        while (std::getline(file, line)) {
            if (line.find("$SOE") != std::string::npos) {
                in_data = true;
                continue;
            }
            if (line.find("$EOE") != std::string::npos) {
                break;
            }
            
            if (in_data && line.find("= A.D.") != std::string::npos) {
                T data;
                if (parseElementSet(file, line, data)) {
                    elements.push_back(data);
                    count++;
                }
            }
        }
        
        std::cout << "Loaded " << count << " " << type << " positions" << std::endl;
        if (!elements.empty()) {
            std::cout << "Time range: JD " << elements.front().jd 
                      << " to " << elements.back().jd << std::endl;
        }
        
        return !elements.empty();
    }
    
    template<typename T>
    bool parseElementSet(std::ifstream& file, const std::string& jd_line, T& data) {
        // Parse JD
        std::istringstream ss(jd_line);
        if (!(ss >> data.jd)) return false;
        
        // Read 4 lines of elements
        std::string line1, line2, line3, line4;
        if (!std::getline(file, line1) || !std::getline(file, line2) || 
            !std::getline(file, line3) || !std::getline(file, line4)) {
            return false;
        }
        
        // Parse elements (simplified - just get what we need)
        if (!extractValue(line1, "EC=", data.eccentricity) ||
            !extractValue(line1, "QR=", data.periapsis_q) ||
            !extractValue(line1, "IN=", data.inclination) ||
            !extractValue(line2, "OM=", data.long_node) ||
            !extractValue(line2, "W =", data.arg_peri) ||
            !extractValue(line3, "TA=", data.true_anomaly)) {
            return false;
        }
        
        // Convert QR from km to AU if needed
        if (data.periapsis_q > 1000.0) {
            data.periapsis_q /= AU;
        }
        
        return true;
    }
    
    bool extractValue(const std::string& line, const std::string& key, double& value) {
        size_t pos = line.find(key);
        if (pos == std::string::npos) return false;
        
        pos += key.length();
        while (pos < line.length() && line[pos] == ' ') pos++;
        
        size_t end = pos;
        while (end < line.length() && line[end] != ' ') end++;
        
        try {
            value = std::stod(line.substr(pos, end - pos));
            return true;
        } catch (...) {
            return false;
        }
    }
    
public:
    bool calculateGeocentricPosition(double target_jd, double& ra, double& dec, double& distance) {
        // Find closest comet data
        const CometData* best = nullptr;
        double min_diff = 1e9;
        
        for (const auto& data : comet_elements) {
            double diff = std::abs(data.jd - target_jd);
            if (diff < min_diff) {
                min_diff = diff;
                best = &data;
            }
        }
        
        if (!best || min_diff > 1.0) {
            std::cout << "No comet data within 1 day of target time" << std::endl;
            return false;
        }
        
        std::cout << "\n=== GEOCENTRIC CALCULATION ===" << std::endl;
        std::cout << "Target JD: " << target_jd << std::endl;
        std::cout << "Using data from JD: " << best->jd << std::endl;
        std::cout << "Time difference: " << (target_jd - best->jd) << " days" << std::endl;
        
        // Step 1: Calculate comet's heliocentric position
        double comet_x, comet_y, comet_z;
        calculateHeliocentricPosition(*best, comet_x, comet_y, comet_z);
        
        std::cout << "\nComet heliocentric position:" << std::endl;
        std::cout << "  [" << comet_x << ", " << comet_y << ", " << comet_z << "] AU" << std::endl;
        
        // Step 2: Calculate Earth's heliocentric position from real data
        double earth_x, earth_y, earth_z;
        if (!calculateEarthPositionFromData(target_jd, earth_x, earth_y, earth_z)) {
            std::cout << "Failed to calculate Earth position from data" << std::endl;
            return false;
        }
        
        std::cout << "\nEarth heliocentric position:" << std::endl;
        std::cout << "  [" << earth_x << ", " << earth_y << ", " << earth_z << "] AU" << std::endl;
        
        // Step 3: Calculate geocentric position (comet - earth)
        double geo_x = comet_x - earth_x;
        double geo_y = comet_y - earth_y;
        double geo_z = comet_z - earth_z;
        
        std::cout << "\nGeocentric position vector:" << std::endl;
        std::cout << "  [" << geo_x << ", " << geo_y << ", " << geo_z << "] AU" << std::endl;
        
        // Step 4: Convert to RA/Dec
        distance = sqrt(geo_x*geo_x + geo_y*geo_y + geo_z*geo_z);
        ra = atan2(geo_y, geo_x) * RAD2DEG;
        dec = asin(geo_z / distance) * RAD2DEG;
        
        // Normalize RA to 0-360
        if (ra < 0) ra += 360.0;
        
        std::cout << "\nFinal geocentric coordinates:" << std::endl;
        std::cout << "  RA: " << std::fixed << std::setprecision(6) << ra << "°" << std::endl;
        std::cout << "  Dec: " << dec << "°" << std::endl;
        std::cout << "  Distance: " << std::setprecision(4) << distance << " AU" << std::endl;
        
        return true;
    }
    
private:
    void calculateHeliocentricPosition(const CometData& data, double& x, double& y, double& z) {
        // Convert orbital elements to Cartesian coordinates
        double nu = data.true_anomaly * DEG2RAD;  // True anomaly
        double r = data.periapsis_q * (1.0 + data.eccentricity * cos(nu)) / 
                   (1.0 + data.eccentricity * cos(nu));  // Distance
        
        // Position in orbital plane
        double x_orb = r * cos(nu);
        double y_orb = r * sin(nu);
        double z_orb = 0.0;
        
        // Rotate to ecliptic coordinates
        double inc = data.inclination * DEG2RAD;
        double node = data.long_node * DEG2RAD;
        double peri = data.arg_peri * DEG2RAD;
        
        // Apply rotation matrices
        double cos_node = cos(node);
        double sin_node = sin(node);
        double cos_inc = cos(inc);
        double sin_inc = sin(inc);
        double cos_peri = cos(peri);
        double sin_peri = sin(peri);
        
        // Combined rotation matrix elements
        double r11 = cos_node * cos_peri - sin_node * sin_peri * cos_inc;
        double r12 = -cos_node * sin_peri - sin_node * cos_peri * cos_inc;
        double r13 = sin_node * sin_inc;
        
        double r21 = sin_node * cos_peri + cos_node * sin_peri * cos_inc;
        double r22 = -sin_node * sin_peri + cos_node * cos_peri * cos_inc;
        double r23 = -cos_node * sin_inc;
        
        double r31 = sin_peri * sin_inc;
        double r32 = cos_peri * sin_inc;
        double r33 = cos_inc;
        
        // Apply rotation
        x = r11 * x_orb + r12 * y_orb + r13 * z_orb;
        y = r21 * x_orb + r22 * y_orb + r23 * z_orb;
        z = r31 * x_orb + r32 * y_orb + r33 * z_orb;
        
        // Convert from ecliptic to equatorial (J2000)
        double obliquity = 23.43929111 * DEG2RAD;
        double x_eq = x;
        double y_eq = y * cos(obliquity) - z * sin(obliquity);
        double z_eq = y * sin(obliquity) + z * cos(obliquity);
        
        x = x_eq;
        y = y_eq;
        z = z_eq;
    }
    
    bool calculateEarthPositionFromData(double target_jd, double& x, double& y, double& z) {
        // Find closest Earth data
        const EarthData* best = nullptr;
        double min_diff = 1e9;
        
        for (const auto& data : earth_elements) {
            double diff = std::abs(data.jd - target_jd);
            if (diff < min_diff) {
                min_diff = diff;
                best = &data;
            }
        }
        
        if (!best || min_diff > 1.0) {
            std::cout << "No Earth data within 1 day of target time" << std::endl;
            return false;
        }
        
        std::cout << "Using Earth data from JD: " << best->jd 
                  << " (diff: " << min_diff << " days)" << std::endl;
        
        // Convert Earth's orbital elements to Cartesian coordinates
        double nu = best->true_anomaly * DEG2RAD;
        double r = best->periapsis_q * (1.0 + best->eccentricity * cos(nu)) / 
                   (1.0 + best->eccentricity * cos(nu));
        
        // Position in orbital plane
        double x_orb = r * cos(nu);
        double y_orb = r * sin(nu);
        double z_orb = 0.0;
        
        // Rotate to ecliptic coordinates
        double inc = best->inclination * DEG2RAD;
        double node = best->long_node * DEG2RAD;
        double peri = best->arg_peri * DEG2RAD;
        
        // Apply rotation matrices (same as for comet)
        double cos_node = cos(node);
        double sin_node = sin(node);
        double cos_inc = cos(inc);
        double sin_inc = sin(inc);
        double cos_peri = cos(peri);
        double sin_peri = sin(peri);
        
        double r11 = cos_node * cos_peri - sin_node * sin_peri * cos_inc;
        double r12 = -cos_node * sin_peri - sin_node * cos_peri * cos_inc;
        double r13 = sin_node * sin_inc;
        
        double r21 = sin_node * cos_peri + cos_node * sin_peri * cos_inc;
        double r22 = -sin_node * sin_peri + cos_node * cos_peri * cos_inc;
        double r23 = -cos_node * sin_inc;
        
        double r31 = sin_peri * sin_inc;
        double r32 = cos_peri * sin_inc;
        double r33 = cos_inc;
        
        x = r11 * x_orb + r12 * y_orb + r13 * z_orb;
        y = r21 * x_orb + r22 * y_orb + r23 * z_orb;
        z = r31 * x_orb + r32 * y_orb + r33 * z_orb;
        
        // Convert from ecliptic to equatorial
        double obliquity = 23.43929111 * DEG2RAD;
        double x_eq = x;
        double y_eq = y * cos(obliquity) - z * sin(obliquity);
        double z_eq = y * sin(obliquity) + z * cos(obliquity);
        
        x = x_eq;
        y = y_eq;
        z = z_eq;
        
        return true;
    }
    
    void calculateEarthPosition(double jd, double& x, double& y, double& z) {
        // Earth's orbital elements (simplified)
        double a = 1.0;  // Semi-major axis (AU)
        double e = 0.0167;  // Eccentricity
        double L0 = 100.464 * DEG2RAD;  // Mean longitude at epoch (J2000)
        
        // Time since J2000
        double T = (jd - 2451545.0) / 365.25;
        
        // Mean longitude
        double L = L0 + 2 * M_PI * T;  // One orbit per year
        
        // Mean anomaly (simplified)
        double M = L - (102.937 * DEG2RAD);  // Subtract longitude of periapsis
        
        // Solve Kepler's equation (first order)
        double E = M + e * sin(M);
        
        // True anomaly
        double nu = 2 * atan(sqrt((1 + e)/(1 - e)) * tan(E/2));
        
        // Distance
        double r = a * (1 - e * cos(E));
        
        // Position in orbital plane
        double x_orb = r * cos(nu);
        double y_orb = r * sin(nu);
        
        // Earth's orbit is nearly in the ecliptic plane
        // Convert to equatorial coordinates
        double obliquity = 23.43929111 * DEG2RAD;
        
        x = x_orb;
        y = y_orb * cos(obliquity);
        z = y_orb * sin(obliquity);
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Precise Geocentric Position Calculator for C/2025 K1" << std::endl;
        std::cout << "Usage: " << argv[0] << " <comet_file> <earth_file> <julian_date>" << std::endl;
        std::cout << "Example: " << argv[0] << " horizons_results-C2025_K1.txt horizons_results-earth.txt 2460770.43" << std::endl;
        return 1;
    }
    
    std::string comet_file = argv[1];
    std::string earth_file = argv[2];
    double target_jd = std::stod(argv[3]);
    
    GeocentricCalculator calc;
    
    std::cout << "=== LOADING HORIZONS DATA ===" << std::endl;
    if (!calc.loadCometFile(comet_file)) {
        std::cout << "Failed to load comet file" << std::endl;
        return 1;
    }
    
    if (!calc.loadEarthFile(earth_file)) {
        std::cout << "Failed to load Earth file" << std::endl;
        return 1;
    }
    
    double ra, dec, distance;
    if (calc.calculateGeocentricPosition(target_jd, ra, dec, distance)) {
        std::cout << "\n=== FINAL RESULT ===" << std::endl;
        std::cout << "Precise geocentric position for JD " << std::fixed << std::setprecision(6) << target_jd << ":" << std::endl;
        std::cout << "RA: " << ra << "°" << std::endl;
        std::cout << "Dec: " << dec << "°" << std::endl;
        std::cout << "Distance: " << std::setprecision(4) << distance << " AU from Earth" << std::endl;
        
        std::cout << "\n=== COMPARISON WITH YOUR FITS ===" << std::endl;
        std::cout << "Your FITS center: RA = 253.829°, Dec = 10.506°" << std::endl;
        std::cout << "Calculated:       RA = " << ra << "°, Dec = " << dec << "°" << std::endl;
        std::cout << "Difference:       RA = " << std::setprecision(3) << (ra - 253.829) 
                  << "°, Dec = " << (dec - 10.506) << "°" << std::endl;
        
        double ra_arcmin = (ra - 253.829) * 60.0;
        double dec_arcmin = (dec - 10.506) * 60.0;
        std::cout << "                  RA = " << ra_arcmin << "', Dec = " << dec_arcmin << "'" << std::endl;
        
        if (std::abs(ra - 253.829) < 1.0 && std::abs(dec - 10.506) < 1.0) {
            std::cout << "\n✓ EXCELLENT MATCH! Precise geocentric calculation working!" << std::endl;
        } else if (std::abs(ra - 253.829) < 5.0 && std::abs(dec - 10.506) < 2.0) {
            std::cout << "\n⚠ Good match - within reasonable precision limits" << std::endl;
        } else {
            std::cout << "\n❌ Large difference - may need time/coordinate system check" << std::endl;
        }
    }
    
    return 0;
}