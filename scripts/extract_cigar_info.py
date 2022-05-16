#!/usr/bin/env python
"""

Given graph alingment with GraphAligner, extract cigar and summarize mutations in the alingment

e.g., 

./extract_cigar_info.py $gaf_file

= match 72950018
X mismatch 202368
I insertions 56026
D deletions 59476




"""

import sys
import re
from collections import defaultdict


mutcol = defaultdict(int)
with open(sys.argv[1]) as infile:
    for line in infile:
        if not line.startswith("#"):
            cigar=line.strip().split()[-1]
            cigar=cigar.split(":")[-1]
            mutsize=""
            mutype=""
            for comp in cigar:
                if re.search("\d+",comp):
                    mutsize = mutsize + comp
                else:
                    mutype=comp
                    mutsize=int(mutsize)
                    mutcol[mutype] += mutsize
                    mutsize=""
                    mutype=""

#CIGAR ID
cigar_id = {

        "M": "matchall",
        "I": "insertions",
        "D": "deletions",
        "N" : "skipped",
        "S" : "softclip",
        "H" : "hardclip",
        "P" : "padding",
        "=" : "match",
        "X" : "mismatch"
        }

for key,value in mutcol.items():
    print(key,cigar_id[key],value)
