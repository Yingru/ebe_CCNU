#!/bin/bash
#script quickcopy.sh
mkdir hydro_stable_2012
cp *.h hydro_stable_2012
cp *.cpp hydro_stable_2012
cp Makefile hydro_stable_2012
cp *.sh hydro_stable_2012
cp *.C hydro_stable_2012
cp *.py hydro_stable_2012
cp *.inp hydro_stable_2012
cp README hydro_stable_2012 
cp Reso3D/ hydro_stable_2012 -r
cp s95p-v1/ hydro_stable_2012 -r
cp SpectraV2/ hydro_stable_2012 -r
cd hydro_stable_2012
mkdir log
mkdir script
mkdir results
mkdir time
