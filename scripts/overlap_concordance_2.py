#!/usr/bin/env python
"""
Will take overlap file from *overlap_coord and calculate concordance from allelic length differ at most x bp 
"""
import sys

foundgeno=0
concor_geno=0
discord_geno=0
missing_1=0
missing_2=0
missing_all=0
allgeno=0

threshold=int(sys.argv[2])

with open(sys.argv[1]) as infile:
    for line in infile:
        token=line.strip().split()
        first_part=token[4:15]
        second_part=token[19:30]
        alldiff=[]
        for al1,al2 in zip(first_part,second_part):
            allgeno += 1
            if al1=="." or al2==".":
                missing_all += 1
                alldiff.append(".")
                if al1==".":
                    missing_1 += 1
                if al2==".":
                    missing_2 += 1
            else:
                alldiff.append(abs(int(al1)-int(al2)))
        #print(*token,*alldiff)

        for comp in alldiff:
            if comp != ".":
                foundgeno += 1
                if comp <= threshold:
                    concor_geno += 1
                if comp > threshold:
                    discord_geno += 1

print(foundgeno,concor_geno,discord_geno,round(concor_geno*100/foundgeno,4),
        missing_1,missing_2,round(missing_1*100/allgeno,4),round(missing_2*100/allgeno,4), 
        missing_all,allgeno,round(missing_all*100/allgeno,4))
