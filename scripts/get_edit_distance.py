#!/usr/bin/env python 
"""

Given gaf alignment, calculate the edit distance from all reads aligned 

"""
import sys

#22_Gaur_6       999999  16304   16797   +       >413    38187   24932   25413   421     496     255     NM:i:75 AS:f:272.5      dv:f:0.15121    id:f:0.84879

contiglen=dict()
totalalign=0
totaledit=0

with open(sys.argv[1]) as infile:
    for line in infile:
        token=line.strip().split()
        contiglen[token[0]]=int(token[1])
        totalalign += abs(int(token[3])-int(token[2]))
        totaledit += int(token[12].split(":")[-1])
        
totalcontig=sum(contiglen.values())
print(totalcontig,totalalign,totaledit)

