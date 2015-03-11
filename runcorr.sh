date
rm -r C12
mkdir C12
g++ Correlation.cpp -o corr
##qsub myjob.sh
ijob=1
while [ "$ijob" -le 8 ]
do
./corr $ijob > corr.out
let "ijob+=1"
done
date
