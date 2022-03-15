#!/usr/bin/env python
"""
Given fasta file, split each entry with size of X bp seq

./split_fa.py input_fa split_size

output: stdout with seqid suffixed with split id

"""
import sys

baselen=int(sys.argv[2])
basecount=0
baseseq=[]
idcount=0

with open(sys.argv[1]) as infile:
    for line in infile:
        if line.startswith(">"):
            suffix=line.strip()[1:]
        else:
            for ind,comp in enumerate(line.strip()):
                if basecount % baselen == 0:
                    if idcount:
                        print(f">{suffix}_{idcount}")
                        print("".join(baseseq))
                    #counter reset
                    idcount += 1
                    basecount += 1
                    baseseq = []
                else:
                    baseseq.append(comp)
                    basecount += 1

    #add the last bases
    print(f">{suffix}_{idcount}")
    print("".join(baseseq))

