#!/bin/bash 
date

ijob=$1
ns=14/$2
while [ "$ijob" -le $2 ]
do
   let "n1=$ns*($ijob-1)+1"
   let "n2=$ns*$ijob"
   echo "#!/bin/bash">>run$ijob.pbs
   echo "">>run$ijob.pbs
   echo "### Job name">>run$ijob.pbs
   echo -n "#PBS -N ppCollision">>run$ijob.pbs
   echo $ijob>>run$ijob.pbs
   echo "#PBS -q lr_batch">>run$ijob.pbs
   echo "#PBS -A ac_jet">>run$ijob.pbs 
   echo "#PBS -l nodes=1:ppn=1">>run$ijob.pbs
   echo "#PBS -l walltime=4:00:00">>run$ijob.pbs
   echo "### Declear job rerunable y/n">>run$ijob.pbs
   echo "#PBS -r n">>run$ijob.pbs
   echo "### Output files">>run$ijob.pbs
   echo -n "#PBS -e log/job">>run$ijob.pbs
   echo -n $ijob>>run$ijob.pbs
   echo ".err">>run$ijob.pbs
   echo -n "#PBS -o log/job">>run$ijob.pbs
   echo -n $ijob>>run$ijob.pbs
   echo ".log">>run$ijob.pbs
   echo "### Mail to user">>run$ijob.pbs
   echo "#PBS -M lgpang.1984@gmail.com">>run$ijob.pbs
   echo "#PBS -m ae">>run$ijob.pbs
   echo "">>run$ijob.pbs
   echo "cd /global/home/users/lpang/TestFreeze/">>run$ijob.pbs
   echo "echo 'Job start at time' ">>run$ijob.pbs
   echo "date">>run$ijob.pbs
   echo "./hydro $n1 $n2 ">>run$ijob.pbs
   echo "echo 'job end at time' ">>run$ijob.pbs
   echo "date">>run$ijob.pbs

   qsub run$ijob.pbs
 
   mv run$ijob.pbs script/.
   let "ijob+=1"
done
