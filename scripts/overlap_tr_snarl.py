#!/usr/bin/env python

import sys

vnid=set()
vnseq=dict()

with open(sys.argv[1]) as infile:
    for line in infile:
        vnid.add(line.strip().split()[3])

with open(sys.argv[2]) as infile:
    for line in infile:
        token=line.strip().split()
        vnseq[f"{token[0]}_{token[1]}_{token[2]}"]=len(token[3])


valid=0
lenlim=0

for input_file in sys.argv[3:]:
    with open(input_file) as infile:
        for line in infile:
            token=line.strip().split()
            if token[0] in vnid:
                missing=sum(len(x.split(",")) <= 2 for x in token[1:])
                if not missing:
                    valid+=1 
                if vnseq[token[0]] >= 10 and vnseq[token[0]] <= 100:
                    lenlim+=1



print(valid)
print(lenlim)
#print(vnid)





