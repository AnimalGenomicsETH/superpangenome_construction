#!/usr/bin/env python

import sys
import re
import math

with open(sys.argv[1]) as infile:
    chromo=sys.argv[1].split("_")[0]
    for line in infile:
        if line[0] != "#":
            svlen=re.findall(r";SVLEN=(-?[0-9]+)",line)
            svoverlap=re.findall(r"SUPP_VEC=([0-1]+)",line)
            print(svlen[0] if svlen else math.nan,svoverlap[0],chromo,sep="\t")


