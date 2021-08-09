#!/usr/bin/env python
"""
This script used to add ref-assembly name in GFA file, it adds "_UCD" in contig name
"""
with open("graph_pggb_chop_test.gfa") as infile, open("gfa_fixed.gfa","w") as outfile:
    for line in infile:
        token=line.strip().split()
        if line.startswith("P"):
            if token[1] in [str(x) for x in range(1,30)]:
                print(token[1])
                token[1]=token[1]+"_UCD"
                print(*token,sep="\t", file=outfile)
            else:
                print(*token,sep="\t", file=outfile) 
        else:
            print(*token,sep="\t", file=outfile)
            