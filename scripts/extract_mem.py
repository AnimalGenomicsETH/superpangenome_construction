#!/usr/bin/env python

import sys
import glob
import re
from collections import defaultdict 


basedir=sys.argv[1]
allfile=glob.glob(basedir+"/*out")
target_stat=["CPU\stime","Max\sMemory","Run\stime"]
resources_id={}
for input_file in allfile:
    chromo=re.findall(r"\d{1,2}",input_file)[0]
    id_res=[]
    with open(input_file) as infile:
        for line in infile:
            if re.search("|".join(target_stat),line):
                #statid=re.findall("|".join(target_stat),line)[0].split()
                id_res.append(float(line.strip().split()[-2]))
                
        if chromo not in resources_id:
            resources_id[chromo]=id_res
        elif resources_id[chromo][0] < id_res[0]:
            resources_id[chromo]=id_res

    id_res=[]

for key,value in resources_id.items():
    print(key,*value,sys.argv[2])




