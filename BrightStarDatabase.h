#include <QDebug>
#include <QString>

class BrightStarDatabase
{
public:
    struct BrightStar {
        QString name;
        QString id;
        double ra, dec, magnitude;
        QString spectralType;
        QString constellation;
    };
    
    static QVector<BrightStar> getAllBrightStars()
    {
        // Comprehensive list of bright stars (magnitude < 3.0) with accurate coordinates
        // Data from the Yale Bright Star Catalog and Hipparcos
        return {
            // Magnitude < 0
            {"Sirius", "HIP_32349", 101.2870, -16.7161, -1.46, "A1V", "Canis Major"},
            {"Canopus", "HIP_30438", 95.9880, -52.6958, -0.74, "A9II", "Carina"},
            {"Arcturus", "HIP_69673", 213.9153, 19.1825, -0.05, "K1.5III", "Boötes"},
            {"Vega", "HIP_91262", 279.2347, 38.7837, 0.03, "A0V", "Lyra"},
            {"Capella", "HIP_24608", 79.1722, 45.9980, 0.08, "G5III", "Auriga"},
            {"Rigel", "HIP_24436", 78.6344, -8.2017, 0.13, "B8Ia", "Orion"},
            {"Procyon", "HIP_37279", 114.8254, 5.2250, 0.34, "F5IV-V", "Canis Minor"},
            {"Achernar", "HIP_7588", 24.6053, -57.2375, 0.46, "B6Vpe", "Eridanus"},
            {"Betelgeuse", "HIP_27989", 88.7929, 7.4071, 0.50, "M1-2Ia", "Orion"},
            {"Hadar", "HIP_68702", 210.9563, -60.3742, 0.61, "B1III", "Centaurus"},
            {"Altair", "HIP_97649", 297.6958, 8.8683, 0.77, "A7V", "Aquila"},
            {"Acrux", "HIP_60718", 186.6497, -63.0992, 0.77, "B0.5IV", "Crux"},
            {"Aldebaran", "HIP_21421", 68.9802, 16.5093, 0.85, "K5III", "Taurus"},
            
            // Magnitude 0-1
            {"Spica", "HIP_65474", 201.2983, -11.1614, 1.04, "B1III-IV", "Virgo"},
            {"Antares", "HIP_80763", 247.3519, -26.4320, 1.09, "M1.5Iab", "Scorpius"},
            {"Pollux", "HIP_37826", 116.3289, 28.0262, 1.14, "K0III", "Gemini"},
            {"Fomalhaut", "HIP_113368", 344.4127, -29.6222, 1.16, "A3V", "Piscis Austrinus"},
            {"Deneb", "HIP_102098", 310.3583, 45.2803, 1.25, "A2Ia", "Cygnus"},
            {"Mimosa", "HIP_62434", 191.9303, -59.6889, 1.30, "B0.5III", "Crux"},
            {"Regulus", "HIP_49669", 152.0930, 11.9672, 1.35, "B7V", "Leo"},
            {"Adhara", "HIP_33579", 104.6564, -28.9720, 1.50, "B2II", "Canis Major"},
            {"Castor", "HIP_36850", 113.6496, 31.8883, 1.57, "A1V", "Gemini"},
            {"Gacrux", "HIP_61084", 187.7915, -57.1133, 1.63, "M3.5III", "Crux"},
            {"Bellatrix", "HIP_25336", 81.2826, 6.3497, 1.64, "B2III", "Orion"},
            {"Elnath", "HIP_25428", 81.5729, 28.6075, 1.68, "B7III", "Taurus"},
            {"Miaplacidus", "HIP_45238", 138.3000, -69.7175, 1.68, "A1III", "Carina"},
            {"Alnilam", "HIP_26311", 84.0533, -1.2019, 1.70, "B0Ia", "Orion"},
            {"Alnair", "HIP_109268", 332.0583, -46.9611, 1.74, "B7IV", "Grus"},
            {"Alioth", "HIP_62956", 193.5073, 55.9598, 1.77, "A1III-IV", "Ursa Major"},
            {"Alnitak", "HIP_26727", 85.1897, -1.9425, 1.79, "O9.7Ibe", "Orion"},
            {"Dubhe", "HIP_54061", 165.9319, 61.7511, 1.79, "K3III", "Ursa Major"},
            {"Mirfak", "HIP_15863", 51.0808, 49.8611, 1.80, "F5Ib", "Perseus"},
            {"Wezen", "HIP_34444", 107.0978, -26.3931, 1.84, "F8Ia", "Canis Major"},
            {"Sargas", "HIP_86228", 264.3297, -42.9986, 1.87, "F1II", "Scorpius"},
            {"Kaus Australis", "HIP_90185", 276.0428, -34.3847, 1.85, "B9.5III", "Sagittarius"},
            {"Avior", "HIP_41037", 125.6283, -59.5097, 1.86, "K3III", "Carina"},
            {"Alkaid", "HIP_67301", 206.8853, 49.3133, 1.86, "B3V", "Ursa Major"},
            
            // Magnitude 1-2
            {"Menkalinan", "HIP_28360", 89.8822, 44.9475, 1.90, "A1IV", "Auriga"},
            {"Atria", "HIP_82273", 253.4153, -69.0278, 1.91, "K2Ib-IIa", "Triangulum Australe"},
            {"Alhena", "HIP_31681", 99.4280, 16.3994, 1.93, "A1.5IV", "Gemini"},
            {"Peacock", "HIP_100751", 306.4119, -56.7350, 1.94, "B2IV", "Pavo"},
            {"Alsephina", "HIP_80816", 247.5553, -25.5919, 2.29, "B2.5V", "Scorpius"},
            {"Mirzam", "HIP_30324", 95.6744, -17.9556, 1.98, "B1II-III", "Canis Major"},
            {"Polaris", "HIP_11767", 37.9544, 89.2642, 1.98, "F7Ib", "Ursa Minor"},
            {"Alphard", "HIP_46390", 141.8967, -8.6586, 1.98, "K3II-III", "Hydra"},
            {"Hamal", "HIP_9884", 31.7933, 23.4625, 2.00, "K2III", "Aries"},
            {"Diphda", "HIP_3419", 10.8975, -17.9867, 2.04, "K0III", "Cetus"},
            {"Nunki", "HIP_92855", 283.8158, -26.2967, 2.02, "B2.5V", "Sagittarius"},
            {"Menkent", "HIP_68933", 211.6708, -36.3700, 2.06, "K0III", "Centaurus"},
            {"Alpheratz", "HIP_677", 2.0970, 29.0906, 2.06, "B8IV", "Andromeda"},
            {"Saiph", "HIP_27366", 86.9392, -9.6697, 2.09, "B0.5Ia", "Orion"},
            {"Algol", "HIP_14576", 47.0421, 40.9564, 2.12, "B8V", "Perseus"},
            {"Almach", "HIP_9640", 30.9747, 42.3297, 2.26, "K3IIb", "Andromeda"},
            {"Eltanin", "HIP_87833", 269.1517, 51.4889, 2.23, "K5III", "Draco"},
            {"Schedar", "HIP_3179", 10.1267, 56.5375, 2.23, "K0IIIa", "Cassiopeia"},
            {"Markab", "HIP_113963", 346.1900, 15.2053, 2.49, "B9III", "Pegasus"},
            {"Denebola", "HIP_57632", 177.2647, 14.5722, 2.14, "A3V", "Leo"},
            {"Navi", "HIP_4427", 14.1772, 60.7175, 2.27, "B0.5IVe", "Cassiopeia"},
            {"Caph", "HIP_746", 2.2950, 59.1500, 2.27, "F2III-IV", "Cassiopeia"},
            {"Izar", "HIP_72105", 221.2467, 27.0744, 2.37, "K0II-III", "Boötes"},
            {"Shaula", "HIP_85927", 263.4022, -37.1039, 1.63, "B1.5IV", "Scorpius"},
            {"Bellatrix", "HIP_25336", 81.2826, 6.3497, 1.64, "B2III", "Orion"},
            {"Naos", "HIP_39429", 120.8958, -40.0033, 2.25, "O4If", "Puppis"},
            {"Mintaka", "HIP_25930", 83.0017, -0.2992, 2.23, "O9.5II", "Orion"},
            {"Arneb", "HIP_25985", 83.1825, -17.8222, 2.58, "F0Ib", "Lepus"},
            {"Gienah", "HIP_102488", 311.1550, 40.2567, 2.46, "K3II", "Cygnus"},
            {"Murzim", "HIP_30324", 95.6744, -17.9556, 1.98, "B1II-III", "Canis Major"},
            {"Algieba", "HIP_50583", 154.9931, 19.8422, 2.28, "K1III", "Leo"},
            {"Merak", "HIP_53910", 165.4600, 56.3825, 2.37, "A1V", "Ursa Major"},
            {"Ankaa", "HIP_2081", 6.5708, -42.3061, 2.39, "K0III", "Phoenix"},
            {"Girtab", "HIP_78401", 239.7119, -19.8050, 2.32, "B0.3IV", "Scorpius"},
            {"Enif", "HIP_107315", 326.0464, 9.8750, 2.39, "K2Ib", "Pegasus"},
            {"Scheat", "HIP_113881", 345.9436, 28.0828, 2.42, "M2.5II-III", "Pegasus"},
            {"Sabik", "HIP_84012", 258.0350, -15.7250, 2.43, "A2V", "Ophiuchus"},
            {"Phecda", "HIP_58001", 178.4567, 53.6947, 2.44, "A0V", "Ursa Major"},
            {"Aludra", "HIP_35904", 111.7897, -29.3033, 2.45, "B5Ia", "Canis Major"},
            {"Markeb", "HIP_45556", 139.2722, -59.2756, 2.47, "B2IV-V", "Vela"},
            {"Nunk", "HIP_92855", 283.8158, -26.2967, 2.02, "B2.5V", "Sagittarius"},
            {"Aspidiske", "HIP_50099", 153.6875, -54.5681, 2.21, "A8Ib", "Carina"}
        };
    }
    
    static QVector<BrightStar> getStarsInField(double centerRA, double centerDec, double radiusDegrees, double maxMagnitude = 3.0)
    {
        QVector<BrightStar> allStars = getAllBrightStars();
        QVector<BrightStar> starsInField;
        
        for (const auto& star : allStars) {
            if (star.magnitude > maxMagnitude) continue;
            
            double ra_diff = star.ra - centerRA;
            double dec_diff = star.dec - centerDec;
            
            // Handle RA wraparound
            if (ra_diff > 180) ra_diff -= 360;
            if (ra_diff < -180) ra_diff += 360;
            
            double cos_dec = cos(centerDec * M_PI / 180.0);
            double angular_distance = sqrt(pow(ra_diff * cos_dec, 2) + pow(dec_diff, 2));
            
            if (angular_distance <= radiusDegrees) {
                starsInField.append(star);
            }
        }
        
        // Sort by magnitude (brightest first)
        std::sort(starsInField.begin(), starsInField.end(), 
                  [](const BrightStar& a, const BrightStar& b) {
                      return a.magnitude < b.magnitude;
                  });
        
        return starsInField;
    }
};
