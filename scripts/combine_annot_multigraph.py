#!/usr/bin/env python
"""

Merged annotation of jasmine SV with the allelic length from graphs 

"""

import sys
import re
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description = __doc__,
                    formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a","--annot",help="annot merged")
    parser.add_argument("-d","--database",help="database as the order of merged",nargs="+")
    return parser.parse_args()


if __name__ == "__main__":

    args=parse_args()
    #varlen database 

    vardict=defaultdict(dict)
    
    for ind,db in enumerate(args.database):
        with open(db) as infile:
            if ind == 0:
                anim_order=next(infile).strip().split()[3:]
                anim_order=[re.findall(r"[A-Za-z]+",x)[0] for x in anim_order]
                anim_id = anim_order
            else:
                anim_id = next(infile).strip().split()[3:]
                anim_id=[re.findall(r"[A-Za-z]+",x)[0] for x in anim_order]
            for line in infile:
                varid=line.strip().split()[2]
                vardict[ind][varid]={x:y for x,y in zip(anim_id,line.strip().split()[3:])}
    
    print("chromo","startpos","mergeid","mergevector","geneid","pathid",*anim_id)

    with open(args.annot) as infile:
        for line in infile:
            tok=line.strip().split()
            if len(tok) == 5:
                varinfo = tok[2].split("_")
                if len(varinfo) ==2 :
                    prog, varid = varinfo
                if len(varinfo) == 3 :
                    prog, varid, dupstat = varinfo
                # if int(prog) in [0,1,2]:
                # somehow there is mixed up in minigraph id
                try:
                   all_allele=vardict[int(prog)][varid]
                   # print(all_allele)
                except:
                    continue
                allele_ordered = [ all_allele[x] for x in anim_order  ]
                print(*tok,varid,*allele_ordered)

                # print(tok,prog,varid)

