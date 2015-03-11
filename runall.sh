
#see myjob.sh

rm result/*
rm result/event* -r
rm log/*

IEVENT=1
echo "Job start at time"
echo `date`
HOME=$(pwd)

### compile initial condition
cd $HOME
rm *.o
cd initial_condition/ || exit
./driver.sh
mv dt_hydro_input.dat ../

###compile hydro evolution subprogram
cd $HOME
rm *.o
make hydro

./hydro $IEVENT

###compile spectra calculation subprogram before resonance decay
cd $HOME
cd SpectraV2/ || exit
rm *.o
g++ calc_spec_before_reso.cpp CSpectra.cpp -O3 -o calc_spec_before_reso

./calc_spec_before_reso $IEVENT

###compile resonance decay subprogram
cd $HOME
cd Reso3D/ || exit
rm *.o
make reso
./reso $IEVENT

###compile spectra calculation subprogram after resonance decay
cd $HOME
cd SpectraV2/ || exit
rm *.o
g++ calc_spec_after_reso.cpp CSpectra.cpp -O3 -o calc_spec_after_reso
./calc_spec_after_reso $IEVENT

echo "Job end at time"
echo `date`

