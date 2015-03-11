#!/bin/bash
#PBS -N Hydrocpp
#PBS -q lr_batch
#PBS -A ac_jet
#PBS -l nodes=1:ppn=1:lr2
#PBS -l walltime=10:00:00
#PBS -e log/job0.err
#PBS -o log/job0.out
#PBS -M lgpang.1984@gmail.com
#PBS -m e
#PBS -m a
#PBS -r n
module load openmpi/1.4.1-gcc
cd /global/home/users/lpang/hydro_ideal_v2/Reso3D
echo "Job start at time"
echo `date`
./reso 5 
echo "Job end at time"
echo `date`
