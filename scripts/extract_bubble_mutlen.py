#!/usr/bin/env python
"""
Take minigraph bubble output and calculate mutation length (alt_len-ref_len)

"""

import sys

with open(sys.argv[1]) as infile:
    for line in infile:
        if line.startswith("chromo"):
            print(line.strip())
        else:
            token=line.strip().split()
            first_part=token[:5]
            second_part=token[6:]
            ref_len=int(token[5].split(":")[1])
            all_len=[]
            missing=0
            for comp in second_part:
                if comp==".":
                    all_len.append(".")
                    missing += 1
                else:
                    alt_len=int(comp.split(":")[1])
                    all_len.append(alt_len-ref_len)
            #print(first_part[0], first_part[1],f"{first_part[3]}{first_part[4]}",*all_len,missing,round(missing/len(second_part),2))
            print(first_part[0], first_part[1],f"{first_part[3]}{first_part[4]}",ref_len,*all_len)
