#!/usr/bin/env python

import sys
import re
from collections import defaultdict


outfile=open(sys.argv[2],"w")
combfile=open(sys.argv[3],"w")

coord_overlap=defaultdict(list)
for db in ["cm","cp","pm","cpm"]:
    with open(f"overlap_{db}_label.tsv") as infile:
        for line in infile:
            token=line.strip().split()
            for comp in token:
                if re.search(r"\d+",comp):
                    coord=int(comp)
                    coord_overlap[coord].append(db)
            

def allele_processor(allele):
    comp=[]
    for alt in allele.split(","):
        if len(alt) >= 50:
            comp.append(alt)
    return(comp)

with open(sys.argv[1]) as infile:
    for line in infile:
        if line.startswith("#"):
            print(line.strip())
        else:
            token=line.strip().split()
            pos=int(token[1])
            token[6]=coord_overlap.get(pos,["private"])

            if token[6] == ["private"]:
                alid=f">{token[1]}"
                for ind,comp in enumerate(allele_processor(token[3]) + allele_processor(token[4])):
                    print(f"{alid}_{ind}",file=outfile)
                    print(comp,file=outfile)

            if "cpm" in coord_overlap.get(pos,["private"]):
                alid=f">{token[1]}"
                for ind,comp in enumerate(allele_processor(token[3]) + allele_processor(token[4])):
                    print(f"{alid}_{ind}",file=combfile)
                    print(comp,file=combfile)

            token[6]=",".join(set(token[6]))
            print(*token)

outfile.close()
combfile.close()
