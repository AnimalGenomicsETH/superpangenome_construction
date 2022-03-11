#!/usr/bin/env python
"""
Add SV support from long read alignment

Usage: ./node_support.py gaf_align_file sv_file > out_file

out_file: similar like sv vcf file but add support (all_allele, ref, min alt, max alt) at the end of each sv

"""

import sys
import re
from collections import defaultdict
from statistics import mean

def node_coverage(input_file):
    nodesup=defaultdict(int)

    with open(input_file) as infile:
        for line in infile:
            allnode=line.strip().split("\t")[5].replace(">","<").split("<")
            for comp in allnode:
                nodesup[comp] += 1
    return nodesup

#coverage based on the edges 
def edges_coverage(input_file):
    edgesup=defaultdict(int)

    with open(input_file) as infile:
        for line in infile:
            allnode=line.strip().split("\t")[5].replace(">","<").split("<")
            for ind,comp in enumerate(allnode[:-1]):
                nodepair=f"{allnode[ind]}_{allnode[ind+1]}"
                edgesup[nodepair] += 1

def allele_checker(allele,nodesup):
    sup = []
    for comp in allele:
        sup.append(nodesup.get(comp,0))
    #return round(mean(sup),2)
    return min(sup)


def edge_checker(allele,edgesup):
    sup = []
    for ind,comp in enumerate(allele[:-1]):
        sup.append(nodesup.get(f"{allele[ind]}_{allele[ind+1]}",0))
    return min(sup)

if __name__ == "__main__":
    nodesup = node_coverage(sys.argv[1])

   #with open(sys.argv[2]) as alignfile:
        #for line in alignfile:
            #if not line.startswith("#"):
                #all_allele = re.findall(r"AT=([^;]+);",line)[0].split(",")
                #allele_sup = []
                #for allele in all_allele:
                #    allnode=allele.replace(">","<").split("<")
                #    allele_sup.append(allele_checker(allnode,nodesup))
                #print(line.strip(),",".join(str(x) for x in allele_sup),allele_sup[0],min(allele_sup[1:]),max(allele_sup[1:]))
                
    with open(sys.argv[2]) as alignfile:
        for line in alignfile:
            if not line.startswith("#"):
                all_allele = re.findall(r"AT=([^;]+);",line)[0].split(",")
                ref_allele = all_allele[0].replace(">","<").split("<")
                ref_sup=allele_checker(ref_allele,nodesup)
                allele_sup = [ref_sup]
                for allele in all_allele[1:]:
                    allnode=allele.replace(">","<").split("<")
                    altnode=[x for x in allnode if x not in ref_allele]
                    if not altnode:
                        allele_sup.append(allele_checker(allnode,nodesup))
                    else:
                        allele_sup.append(allele_checker(altnode,nodesup))
                print(line.strip(),",".join(str(x) for x in allele_sup),allele_sup[0],min(allele_sup[1:]),max(allele_sup[1:]))
