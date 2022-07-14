#!/usr/bin/env python

import sys

with open(sys.argv[1]) as infile:
    anim=next(infile).strip().split()[4:]
    print("UCD",*anim)
    geno_tidy=[0]
    for line in infile:
        for geno in line.strip().split()[9:]:
            if "/" in geno:
                geno_tidy.append(geno.split("/")[0])
            elif "." in geno:
                geno_tidy.append("NA")
            else:
                geno_tidy.append(geno)
        print(*geno_tidy)
        geno_tidy=[0]
