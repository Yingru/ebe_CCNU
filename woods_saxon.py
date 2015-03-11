#!/usr/bin/python
#filename:woods_saxon.py
from math import *
from Numeric import *
import random
import array

def ws(r):
    return 1.0/(exp((r-7.0)/0.54)+1.0)

fmax=ws(0.0)

i=0
ary=zeros(20)

while i<=1000000:
    r=20.0*random.random()
    f=random.random()
    if f<ws(r)/fmax:
       ary[int(r)] += 1 
    i += 1   

for i in range(0,20):
    print i,ary[i]
