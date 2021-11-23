#!/usr/bin/env python

import sys
from collections import defaultdict

dict_count=defaultdict(int)

with open(sys.argv[1]) as infile:
    for line in infile:
        all_anim=line.strip().split()[2].split(",")
        for comp in all_anim:
            dict_count[comp] += 1

print(dict_count)
            


