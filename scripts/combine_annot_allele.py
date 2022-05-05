#!/usr/bin/env python
"""

Merged annotation of jasmine SV with the allelic length from graphs 

"""

import sys
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
            for line in infile:
                varid=line.strip().split()[2]
                vardict[ind][varid]=line.strip().split()[3:]

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
                   print(*tok,varid,*vardict[int(prog)][varid])
                except:
                    continue

                # print(tok,prog,varid)

