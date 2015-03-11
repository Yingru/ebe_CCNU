rm *.o
g++ calc_spec_before_reso.cpp CSpectra.cpp -O3 -o calc_spec_before_reso

rm *.o
g++ calc_spec_after_reso.cpp CSpectra.cpp -O3 -o calc_spec_after_reso

