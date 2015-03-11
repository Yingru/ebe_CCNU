#!/bin/bash

ijob=$1
while [ "$ijob" -le $2 ]
do 
qdel $ijob
let "ijob+=1"
done
