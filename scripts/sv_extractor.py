#!/usr/bin/env python

import sys

def allele_processor(allele):
    comp=[]
    for alt in allele.split(","):
        if len(alt) >= 50:
            comp.append(alt)
    return(comp)



with open(sys.argv[1]) as infile:
    for line in infile:
        if not line.startswith("#"):
            token=line.strip().split()
            pos=int(token[1])
            alid=f">{token[1]}"
            for ind,comp in enumerate(allele_processor(token[3]) + allele_processor(token[4])):
                 print(f"{alid}_{ind}")
                 print(comp)



