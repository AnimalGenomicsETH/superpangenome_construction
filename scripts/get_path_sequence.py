#!/usr/bin/env python
import sys
import os 
import re
from collections import defaultdict

def get_nonref(graphfile,refnode):
    nodedict=dict()
    nodeseq=dict()
    totalref=0
    nonreflen=defaultdict(int)
    with open(graphfile) as infile:
        for line in infile:
            token=line.strip().split()
            if line[0] == "S":
                nodedict[token[1]] = len(token[2])
                nodeseq[token[1]] = token[2]
            elif line[0] == "P" and re.search("UCD",token[1]):
                nodeall=token[2].split(",")
                for node in nodeall:
                    totalref += nodedict[node[:-1]]
            elif line[0] == "P" and not re.search("UCD",token[1]):
                nodeall=token[2].split(",")
                anim=re.findall(r'[A-Za-z]+',token[1])[0]
                seqnonref=""
                seqid=1
                prevnonref=0
                for node in nodeall:
                    node=node[:-1]
                    if node not in refnode:
                        nonreflen[anim] += nodedict[node]
                        if prevnonref:
                            seqnonref=seqnonref + nodeseq[node]
                            prevnonref=1
                        else:
                            if len(seqnonref) >= 50:
                                print(f">{anim}_{seqid}")
                                print(seqnonref)
                            seqnonref=nodeseq[node]
                            seqid+=1
                            prevnonref=0
    return [totalref, nonreflen]


def get_refnode(input_file):
    with open(input_file) as infile:
        for line in infile:
            if line[0] == "P":
                if re.search("UCD",line):
                    node_all = [x[:-1] for x in line.strip().split()[2].split(",")]
    return set(node_all)


if __name__ == "__main__":
    refnode=get_refnode(sys.argv[1])
    trs=get_nonref(sys.argv[1],refnode)






