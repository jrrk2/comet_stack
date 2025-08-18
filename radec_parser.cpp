#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>

struct ObservationData {
    double jd;
    double ra;   // degrees
    double dec;  // degrees
    std::string visibility; // visibility flags (m, *, C, N, A, etc.)
};

class RADecParser {
private:
    std::vector<ObservationData> observations;
    
public:
    bool loadFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "Cannot open file: " << filename << std::endl;
            return false;
        }
        
        std::cout << "Loading RA/Dec file: " << filename << std::endl;
        
        std::string line;
        bool in_data = false;
        int count = 0;
        
        while (std::getline(file, line)) {
            if (line.find("$$SOE") != std::string::npos) {
                in_data = true;
                continue;
            }
            if (line.find("$$EOE") != std::string::npos) {
                break;
            }
            
            if (in_data && !line.empty() && line[0] == ' ') {
                ObservationData obs;
                if (parseLine(line, obs)) {
                    observations.push_back(obs);
                    count++;
                }
            }
        }
        
        std::cout << "Loaded " << count << " observations" << std::endl;
        if (!observations.empty()) {
            std::cout << "Time range: JD " << std::fixed << std::setprecision(3)
                      << observations.front().jd << " to " << observations.back().jd << std::endl;
        }
        
        return !observations.empty();
    }
    
private:
    bool parseLine(const std::string& line, ObservationData& obs) {
        // Format: " 2024-Jan-01 00:00 2460310.500000000  m  324.99787  14.79833"
        std::istringstream ss(line);
        
        std::string date, time, jd_str, vis_flag, ra_str, dec_str;
        
        if (!(ss >> date >> time >> jd_str)) {
            return false;
        }
        
        // Parse JD
        try {
            obs.jd = std::stod(jd_str);
        } catch (...) {
            return false;
        }
        
        // Get visibility flag (optional)
        if (ss >> vis_flag) {
            obs.visibility = vis_flag;
            
            // Parse RA and Dec
            if (ss >> ra_str >> dec_str) {
                try {
                    obs.ra = std::stod(ra_str);
                    obs.dec = std::stod(dec_str);
                    return true;
                } catch (...) {
                    return false;
                }
            }
        }
        
        return false;
    }
    
public:
    bool getPosition(double target_jd, double& ra, double& dec) {
        if (observations.empty()) {
            std::cout << "No observations loaded" << std::endl;
            return false;
        }
        
        // Find the closest observation by time
        auto closest = std::min_element(observations.begin(), observations.end(),
            [target_jd](const ObservationData& a, const ObservationData& b) {
                return std::abs(a.jd - target_jd) < std::abs(b.jd - target_jd);
            });
        
        if (closest == observations.end()) {
            return false;
        }
        
        double time_diff = target_jd - closest->jd;
        
        std::cout << "\n=== POSITION LOOKUP ===" << std::endl;
        std::cout << "Target JD: " << std::fixed << std::setprecision(6) << target_jd << std::endl;
        std::cout << "Found data at JD: " << closest->jd << std::endl;
        std::cout << "Time difference: " << std::setprecision(4) << time_diff << " days" << std::endl;
        std::cout << "Time difference: " << std::setprecision(1) << (time_diff * 24.0) << " hours" << std::endl;
        
        if (std::abs(time_diff) > 0.1) {
            std::cout << "⚠ Warning: Large time difference - consider interpolation" << std::endl;
            
            // Try simple linear interpolation if we have nearby points
            if (tryInterpolation(target_jd, ra, dec)) {
                return true;
            }
        }
        
        ra = closest->ra;
        dec = closest->dec;
        
        std::cout << "Visibility: " << closest->visibility << std::endl;
        std::cout << "Position: RA = " << ra << "°, Dec = " << dec << "°" << std::endl;
        
        return true;
    }
    
private:
    bool tryInterpolation(double target_jd, double& ra, double& dec) {
        // Find the two closest points that bracket the target time
        ObservationData* before = nullptr;
        ObservationData* after = nullptr;
        
        for (auto& obs : observations) {
            if (obs.jd <= target_jd) {
                if (!before || obs.jd > before->jd) {
                    before = &obs;
                }
            }
            if (obs.jd >= target_jd) {
                if (!after || obs.jd < after->jd) {
                    after = &obs;
                }
            }
        }
        
        if (before && after && before != after) {
            double dt = after->jd - before->jd;
            double t = (target_jd - before->jd) / dt;
            
            // Linear interpolation
            ra = before->ra + t * (after->ra - before->ra);
            dec = before->dec + t * (after->dec - before->dec);
            
            // Handle RA wraparound at 0/360
            double ra_diff = after->ra - before->ra;
            if (ra_diff > 180.0) {
                ra = before->ra + t * (ra_diff - 360.0);
            } else if (ra_diff < -180.0) {
                ra = before->ra + t * (ra_diff + 360.0);
            }
            
            // Normalize RA
            if (ra < 0) ra += 360.0;
            if (ra >= 360.0) ra -= 360.0;
            
            std::cout << "Interpolated between JD " << before->jd << " and " << after->jd << std::endl;
            std::cout << "Interpolation factor: " << std::setprecision(3) << t << std::endl;
            
            return true;
        }
        
        return false;
    }
    
public:
    void showTimeRange() {
        if (observations.empty()) return;
        
        std::cout << "\n=== TIME COVERAGE ===" << std::endl;
        std::cout << "Total observations: " << observations.size() << std::endl;
        
        // Check coverage around August 16, 2025
        double target_jd = 2460770.43; // Your observation
        
        auto closest = std::min_element(observations.begin(), observations.end(),
            [target_jd](const ObservationData& a, const ObservationData& b) {
                return std::abs(a.jd - target_jd) < std::abs(b.jd - target_jd);
            });
        
        if (closest != observations.end()) {
            double diff = std::abs(closest->jd - target_jd);
            std::cout << "Closest to your observation: JD " << closest->jd 
                      << " (diff: " << std::setprecision(4) << diff << " days)" << std::endl;
            
            if (diff < 0.5) {
                std::cout << "✓ Excellent coverage for your observation!" << std::endl;
            } else if (diff < 5.0) {
                std::cout << "⚠ Reasonable coverage, interpolation recommended" << std::endl;
            } else {
                std::cout << "❌ Poor coverage for your observation date" << std::endl;
            }
        }
        
        // Show first and last few entries
        std::cout << "\nFirst few observations:" << std::endl;
        for (size_t i = 0; i < std::min(size_t(3), observations.size()); i++) {
            const auto& obs = observations[i];
            std::cout << "  JD " << obs.jd << ": RA=" << obs.ra << "°, Dec=" << obs.dec << "°" << std::endl;
        }
        
        std::cout << "\nLast few observations:" << std::endl;
        size_t start = std::max(size_t(0), observations.size() - 3);
        for (size_t i = start; i < observations.size(); i++) {
            const auto& obs = observations[i];
            std::cout << "  JD " << obs.jd << ": RA=" << obs.ra << "°, Dec=" << obs.dec << "°" << std::endl;
        }
    }
    
    void generateEphemeris(double start_jd, double end_jd, double step_hours, 
                          const std::string& output_file) {
        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cout << "Cannot create output file: " << output_file << std::endl;
            return;
        }
        
        out << "# C/2025 K1 (ATLAS) Geocentric Ephemeris" << std::endl;
        out << "# Observer location: 0.17°E, 52.19°N (your coordinates)" << std::endl;
        out << "JD,RA(deg),Dec(deg),RA(hms),Dec(dms)" << std::endl;
        
        double step_days = step_hours / 24.0;
        int count = 0;
        
        for (double jd = start_jd; jd <= end_jd; jd += step_days) {
            double ra, dec;
            if (getPosition(jd, ra, dec)) {
                // Convert RA to hours:minutes:seconds
                double ra_hours = ra / 15.0;
                int h = static_cast<int>(ra_hours);
                int m = static_cast<int>((ra_hours - h) * 60);
                double s = ((ra_hours - h) * 60 - m) * 60;
                
                // Convert Dec to degrees:minutes:seconds
                int d = static_cast<int>(dec);
                int am = static_cast<int>(std::abs((dec - d) * 60));
                double as = (std::abs((dec - d) * 60) - am) * 60;
                char dec_sign = (dec >= 0) ? '+' : '-';
                
                out << std::fixed << std::setprecision(6) << jd << ","
                    << ra << "," << dec << ","
                    << h << ":" << std::setfill('0') << std::setw(2) << m << ":" 
                    << std::setprecision(2) << s << ","
                    << dec_sign << std::setfill('0') << std::setw(2) << std::abs(d) << ":" 
                    << std::setw(2) << am << ":" << std::setprecision(1) << as << std::endl;
                count++;
            }
        }
        
        std::cout << "Generated " << count << " ephemeris points in " << output_file << std::endl;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "RA/Dec Horizons Parser for C/2025 K1" << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "  " << argv[0] << " info <radec_file>" << std::endl;
        std::cout << "  " << argv[0] << " pos <radec_file> <julian_date>" << std::endl;
        std::cout << "  " << argv[0] << " ephemeris <radec_file> <start_jd> <end_jd> <step_hours> <output.csv>" << std::endl;
        std::cout << std::endl;
        std::cout << "Examples:" << std::endl;
        std::cout << "  " << argv[0] << " info ~/Downloads/horizons_results-C2025_K1-radec.txt" << std::endl;
        std::cout << "  " << argv[0] << " pos ~/Downloads/horizons_results-C2025_K1-radec.txt 2460770.43" << std::endl;
        std::cout << "  " << argv[0] << " ephemeris ~/Downloads/horizons_results-C2025_K1-radec.txt 2460770 2460771 1 ephemeris.csv" << std::endl;
        return 1;
    }
    
    std::string command = argv[1];
    std::string radec_file = argv[2];
    
    RADecParser parser;
    
    if (command == "info") {
        if (parser.loadFile(radec_file)) {
            parser.showTimeRange();
        }
        
    } else if (command == "pos" && argc >= 4) {
        double jd = std::stod(argv[3]);
        
        if (parser.loadFile(radec_file)) {
            double ra, dec;
            if (parser.getPosition(jd, ra, dec)) {
                std::cout << "\n=== RESULT ===" << std::endl;
                std::cout << "Geocentric position for JD " << std::fixed << std::setprecision(6) << jd << ":" << std::endl;
                std::cout << "RA: " << ra << "°" << std::endl;
                std::cout << "Dec: " << dec << "°" << std::endl;
                
                std::cout << "\n=== COMPARISON WITH YOUR FITS ===" << std::endl;
                std::cout << "Your FITS center: RA = 253.829°, Dec = 10.506°" << std::endl;
                std::cout << "Horizons geocentric: RA = " << ra << "°, Dec = " << dec << "°" << std::endl;
                std::cout << "Difference: RA = " << std::setprecision(3) << (ra - 253.829) 
                          << "°, Dec = " << (dec - 10.506) << "°" << std::endl;
                
                double ra_arcmin = (ra - 253.829) * 60.0;
                double dec_arcmin = (dec - 10.506) * 60.0;
                std::cout << "            RA = " << ra_arcmin << "', Dec = " << dec_arcmin << "'" << std::endl;
                
                if (std::abs(ra - 253.829) < 0.5 && std::abs(dec - 10.506) < 0.5) {
                    std::cout << "\n✓ EXCELLENT MATCH! Horizons data matches your FITS perfectly!" << std::endl;
                } else if (std::abs(ra - 253.829) < 2.0 && std::abs(dec - 10.506) < 1.0) {
                    std::cout << "\n⚠ Good match - small differences likely due to time/coordinate precision" << std::endl;
                } else {
                    std::cout << "\n❌ Significant difference - check observation time or coordinate systems" << std::endl;
                }
            }
        }
        
    } else if (command == "ephemeris" && argc >= 7) {
        double start_jd = std::stod(argv[3]);
        double end_jd = std::stod(argv[4]);
        double step_hours = std::stod(argv[5]);
        std::string output_file = argv[6];
        
        if (parser.loadFile(radec_file)) {
            parser.generateEphemeris(start_jd, end_jd, step_hours, output_file);
        }
        
    } else {
        std::cout << "Invalid command or missing arguments" << std::endl;
        return 1;
    }
    
    return 0;
}