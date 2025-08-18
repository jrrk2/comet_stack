g++ ../radec_parser.cpp -o radec_parser
./radec_parser ephemeris ../horizons_C2025_K1.txt 2460903.8 2460904.1 0.5 comet_ephemeris.csv\n
cmake ..
make
./comet_stacker comet_ephemeris.csv ~/Downloads/C2025_K1_\(ATLAS\)_2025-08-16_22-24-37 comet_stack.fits
