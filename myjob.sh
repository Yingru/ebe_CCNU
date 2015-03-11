#!/bin/bash
#PBS -N hydro3D
#PBS -q lr_batch
#PBS -A ac_jet
#PBS -l nodes=1:ppn=1:lr2
#PBS -l walltime=10:00:00
#PBS -e log/job0.err
#PBS -o log/job0.out
#PBS -M lpang@lbl.gov
#PBS -m e
#PBS -m a
#PBS -r n

#module load openmpi/1.4.1-gcc

rm log/*

IEVENT=1
echo "Job start at time"
echo `date`
###compile hydro evolution subprogram
HOME=$(pwd)

#cd $HOME
#rm *.o
#make hydro
#
#./hydro $IEVENT

###compile spectra calculation subprogram before resonance decay
cd $HOME
cd SpectraV2/
rm *.o
g++ calc_spec_before_reso.cpp CSpectra.cpp -O3 -o calc_spec_before_reso

./calc_spec_before_reso $IEVENT

###compile resonance decay subprogram
cd $HOME
cd Reso3D/
rm *.o
make reso
./reso $IEVENT

###compile spectra calculation subprogram after resonance decay
cd $HOME
cd SpectraV2/
rm *.o
g++ calc_spec_after_reso.cpp CSpectra.cpp -O3 -o calc_spec_after_reso
./calc_spec_after_reso $IEVENT

echo "Job end at time"
echo `date`
