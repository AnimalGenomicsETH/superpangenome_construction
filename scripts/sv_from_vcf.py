#!/usr/bin/env python
    
import sys

with open(sys.argv[1]) as infile:
    for line in infile:
        if line[0] == "#":
            print(line.strip())
        else:
            token=line.strip().split()
            if len(token[3]) >= 50:
                print(line.strip())
                continue
            for alt in token[4].split(","):
                if len(alt) >= 50:
                    print(line.strip())
                    continue

            
            
