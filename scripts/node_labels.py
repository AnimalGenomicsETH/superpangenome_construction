#!/usr/bin/env python

import argparse
import re
from dataclasses import dataclass, field
from typing import List

def parse_args():
    parser = argparse.ArgumentParser(description = __doc__,
                    formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--input",help="input graph file")
    parser.add_argument("-o","--output",help="output node label file")
    return parser.parse_args()

if __name__ == "__main__":
    args=parse_args()
    input=args.input
    output=args.output
    
    @dataclass
    class node_info:
        name:str
        length:int=0
        label:List[str] = field(default_factory=list)
    
    node_comb=dict()
    
    with open(input) as infile:
        for line in infile:
            token=line.strip().split()
            if token[0] == "S":
                node_id=token[1]
                node_len=len(token[2])
                node_comb[node_id]=node_info(name=node_id,length=node_len)
            elif token[0] == "P":
                label_id = re.findall(pattern=r'[A-Za-z]+',string=token[1])[0]
                label_node = token[2].split(",")
                for comp in label_node:
                    comp = comp[:-1]
                    node_comb[comp].label.append(label_id)

    with open(output,"w") as outfile:
        outfile.write("node_id\tnode_len\tnode_label\tshort_label\n")
        for key,value in node_comb.items():
            short_label="".join(x[0] for x in sorted(value.label)) 
            print(value.name,value.length,",".join(sorted(value.label)),short_label,sep="\t",file=outfile)



    

    
