#!/usr/bin/env python 
"""

Given gaf alignment, calculate the edit distance from all reads aligned 

"""
import sys
import re
import os
from collections import defaultdict

#22_Gaur_6       999999  16304   16797   +       >413    38187   24932   25413   421     496     255     NM:i:75 AS:f:272.5      dv:f:0.15121    id:f:0.84879

contiglen=dict()
totalalign=0
totaledit=0
binstat=defaultdict(lambda: [0,0,0,0])

with open(sys.argv[1]) as infile:
    for line in infile:
        token=line.strip().split()
        contiglen[token[0]]=int(token[1])
        totalalign += abs(int(token[3])-int(token[2]))
        totaledit += int(token[12].split(":")[-1])
        #this is per bin stat
        chromo=token[0].split("_")[0]
        contigid=token[0].split("_")[-1]
        binstat[contigid][0]=chromo
        binstat[contigid][1]=int(token[1])
        binstat[contigid][2]+=abs(int(token[3])-int(token[2]))
        binstat[contigid][3]+=int(token[12].split(":")[-1])


## id of the file
##26_cactus_BS8.gaf
chromo,prog,animid=os.path.basename(sys.argv[1]).split(".")[0].split("_")

for key,value in binstat.items():
    print(chromo,prog,animid,key,*value)


        
# totalcontig=sum(contiglen.values())
# print(totalcontig,totalalign,totaledit)

