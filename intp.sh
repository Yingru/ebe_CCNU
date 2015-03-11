#!/bin/bash
#PBS -N intp
#PBS -q lr_batch
#PBS -A ac_jet
#PBS -l nodes=1:ppn=1:lr2
#PBS -l walltime=00:10:00
#PBS -e log/job0.err
#PBS -o log/job0.out
#PBS -M lpang@lbl.gov
#PBS -m e
#PBS -m a
#PBS -r n
cd /global/home/users/lpang/hydro_b0-10
g++ Interpolation.cpp
echo "Job start at time"
echo `date`
./a.out
echo "Job end at time"
echo `date`
