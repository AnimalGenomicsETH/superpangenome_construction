#!/usr/bin/env python
"""
Graph minigraph call file, output the non-ref sequence for particular breed

"""
import sys 
import re

def get_node_seq(graphfile):
    nodeseq=dict()
    with open(graphfile) as infile:
        for line in infile:
            if line[0] == "S":
                token=line.strip().split()
                nodeseq[token[1]] = token[2]
    return nodeseq


if __name__ == "__main__":
    nodeseq=get_node_seq(sys.argv[1])

    with open(sys.argv[2]) as infile:
        seen_node=set()
        anim=""
        for line in infile:
            token=line.strip().split()
            chromo,start,stop,startnode,stopnode,bubble=token
            if "." not in bubble and "*" not in bubble:
                bubble_comp=bubble.split(":")
                if not anim:
                    anim=re.findall(pattern=r"[a-zA-Z]+",string=bubble)[-1]
                bubble_node=bubble_comp[0].replace("<",">").split(">")
                bubble_seq=""
                for comp in bubble_node:
                    if comp:
                        bubble_seq += nodeseq[comp]
                if len(bubble_seq) >= 50:
                    print(f">{anim}_{chromo}_{start}_{stop}")
                    print(bubble_seq)
            bubble_seq=""



