#!/usr/bin/env python
"""
This script will take the minigraph combined bubble bed
and output the label for each node

"""
import sys
from collections import defaultdict

node_label=defaultdict(list)

with open(sys.argv[1]) as infile:
    next(infile)
    for line in infile:
        token=line.strip().split()[5:]
        for comp in token:
            if comp != ".":
                node_list=comp.split(":")[0].replace("<",">").split(">")
                breed=comp.split(":")[3].split("_")[1]
                for node in node_list:
                    node_label[node].append(breed)

print("node_ID\tlong_label\tref_status\tshort_label\tnode_stat")
for key,values in node_label.items():
    if key:
        if key != "*":
            long_label=",".join(values)
            ref_status = "R" if "UCD" in values else "NR"
            short_label="".join(x[0] for x in sorted(values))
            node_stat="R" if "UCD" in values else short_label
            print(key,long_label,ref_status,short_label,node_stat,sep="\t")

