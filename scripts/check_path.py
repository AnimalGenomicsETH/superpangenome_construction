#!/usr/bin/env python
"""

Given assembly alignment to graph, check whether the GAF alignment conforms with the path in the graph 

./check_path.py -g {graph file} -a {GAF assembly alignment to matched graph} > outfile 

stdout will output:

path edge_alignment_inpath edge_alignment_outpath nodes_notpart_of_path

"""

import argparse
import re
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description = __doc__,
                    formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g","--graph",help="input graph file")
    parser.add_argument("-a","--align",help="alignment file")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    graph = args.graph
    align = args.align 
    edge_comb = defaultdict(lambda: defaultdict(set))
    with open(graph) as infile:
        for line in infile:
            if line.startswith("P"):
                token=line.strip().split()
                breed=re.findall(r"[A-Za-z]+",token[1])[0]
                nodes_list = token[2].split(",")
                for ind,comp in enumerate(nodes_list[:-1]):
                    #parent, child = sorted([int(comp[:-1]), int(nodes_list[ind+1][:-1])])
                    parent, child = [int(comp[:-1]), int(nodes_list[ind+1][:-1])]
                    edge_comb[breed][parent].add(child)
                    edge_comb[breed][child].add(parent)
    for key, value in edge_comb.items():
        for comb1, comb2 in value.items():
            print(key,comb1,comb2)

    #in_path, out_path, not_found 
    support_res = defaultdict(lambda: [0,0,0])
    with open(align) as infile:
        for line in infile:
            token=line.strip().split()
            breed = re.findall(r"[A-Za-z]+",token[0])[0]
            node_comp = token[5].replace(">","<").split("<")
            for ind,comp in enumerate(node_comp[:-1]):
                if comp:
                    parent = int(comp) 
                    child = int(node_comp[ind+1])
                    print(parent,child)
                    if edge_comb[breed].get(parent,0) and edge_comb[breed].get(child,0):
                        if edge_comb[parent].get(child,0):
                            support_res[breed][0] += 1
                        elif edge_comb[child].get(parent,0):
                            support_res[breed][0] += 1
                        else:
                            support_res[breed][1] += 1
                    else:
                        support_res[breed][2] += 1
    
   
    for key, value in support_res.items():
        print(key,*value)

