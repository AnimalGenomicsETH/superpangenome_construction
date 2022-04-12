#!/usr/bin/env python
"""

Given the snarl deconstrution vcf, output the frequency of alelle (paths)
- Not including missing 
- Not including multiple paths per assembly, if it is not unique

"""

import sys
import re
from collections import defaultdict

with open(sys.argv[1]) as infile:
    for line in infile:
        if line[0] != "#":
            token=line.strip().split()
            geno=token[9:]
            altfreq=defaultdict(int)
            missing=0
            total_anim=len(geno)
            total_allele=total_anim*2
            for comp in geno:
                if re.search(r"\.",comp):
                    missing += 1
                elif re.search("/",comp):
                    all_geno = list(set(comp.split("/")))
                    if len(all_geno)==1 and all_geno[0] != "0":
                        altfreq[all_geno[0]] += 1
                else:
                    if comp != "0":
                        altfreq[comp] += 1
            for key,value in altfreq.items():
                print(token[0],token[1],key,value,missing,total_anim,total_allele,(total_anim-missing)*2)
                        


