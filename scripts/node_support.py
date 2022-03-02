#!/usr/bin/env python
"""
Add SV support from long read alignment

Usage: ./node_support.py gaf_align_file sv_file > out_file

out_file: similar like sv vcf file but add support (all_allele, ref, min alt, max alt) at the end of each sv

"""

import sys
import re
from collections import defaultdict


def node_coverage(input_file):
    nodesup=defaultdict(int)

    with open(input_file) as infile:
        for line in infile:
            allnode=line.strip().split("\t")[5].replace(">","<").split("<")
            for comp in allnode:
                nodesup[comp] += 1
    return nodesup


def allele_checker(allele,nodesup):
    sup = []
    for comp in allele:
        sup.append(nodesup.get(comp,0))
    return min(sup)


if __name__ == "__main__":
    nodesup = node_coverage(sys.argv[1])

    with open(sys.argv[2]) as alignfile:
        for line in alignfile:
            if not line.startswith("#"):
                all_allele = re.findall(r"AT=([^;]+);",line)[0].split(",")
                allele_sup = []
                for allele in all_allele:
                    allnode=allele.replace(">","<").split("<")
                    allele_sup.append(allele_checker(allnode,nodesup))
                print(line.strip(),",".join(str(x) for x in allele_sup),allele_sup[0],min(allele_sup[1:]),max(allele_sup[1:]))
