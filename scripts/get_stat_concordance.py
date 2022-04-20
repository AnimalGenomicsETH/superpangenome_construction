#!/usr/bin/env python
"""
Calculate concordance stat between query and truth

Usage:
    ./get_stat_concordance.py query.vcf(gzip) truth.vcf(gzip)

"""

import sys
import re
import gzip
from collections import defaultdict


#if more than two files: Combining concordance stat from single chromosome results
if len(sys.argv) > 3:
    concor_table = defaultdict(lambda: [0 for x in range(0,12)])
    for input_file in sys.argv[1:]:
        with open(input_file) as infile:
            for line in infile:
                token=line.strip().split()
                if len(token) > 10:
                    anim_id=token[0]
                    for ind,comp in enumerate(token[1:]):
                        concor_table[anim_id][ind] += int(comp)
    # Concordance calculation
    all_geno_match = concor_table["All"]

    # concordance stat
    print("Concordance", sum(all_geno_match[x] for x in [0,4,8])/sum(all_geno_match[x] for x in range(0,9)))
    print("NR Sensitivity", sum(all_geno_match[x] for x in [4,5,7,8])/sum(all_geno_match[x] for x in [1,2,4,5,7,8,10,11]))
    print("NR Discrepancy", sum(all_geno_match[x] for x in [2,3,5,6,7])/sum(all_geno_match[x] for x in range(1,9)))
    print("Precision", sum(all_geno_match[x] for x in [4,8])/sum(all_geno_match[x] for x in range(3,9)))

    sys.exit(0)

query=gzip.open(sys.argv[1],"rt")
truth=gzip.open(sys.argv[2],"rt")

#get the id of anim 

def get_id(infile):
    for line in infile:
        if line[:6] == "#CHROM":
            return(line.strip().split("\t")[9:])

query_id=get_id(query)
truth_id=get_id(truth)

query_line=next(query).strip().split()
truth_line=next(truth).strip().split()

query_pos=query_line[1]
truth_pos=truth_line[1]

# print(query_line[1])
# print(truth_line[1])

def geno_encoder(geno):
    """
    Allelic dosage encoder with possibility of phasing
    """

    if re.search(r"0(\||/)0",geno):
        return 0 
    elif re.search(r"0(\||/)1",geno) or re.search(r"1(\||/)0",geno):
        return 1
    elif re.search(r"1(\||/)1",geno):
        return 2
    else:
        return 3
     

def geno_dict(vcfline,anim_id):
    """
    return dict of anim_id: allelic_dosage
    """

    return({ x:geno_encoder(y.split(":")[0]) for x,y in zip(anim_id,vcfline[9:]) })

def concor_match(query_geno,truth_geno):
    """

    Index of match query and truth as in the paper

        truth
           0 1  2 
         0 0 1  2 
         1 3 4  5
    test 2 6 7  8   
         3 9 10 11
         
    """
    
    pair_geno = [query_geno,truth_geno]
    if pair_geno == [0,0]:
        return 0
    elif pair_geno == [0,1]:
        return 1
    elif pair_geno == [0,2]:
        return 2
    elif pair_geno == [1,0]:
        return 3
    elif pair_geno == [1,1]:
        return 4
    elif pair_geno == [1,2]:
        return 5
    elif pair_geno == [2,0]:
        return 6
    elif pair_geno == [2,1]:
        return 7
    elif pair_geno == [2,2]:
        return 8
    elif pair_geno == [3,0]:
        return 9
    elif pair_geno == [3,1]:
        return 10
    elif pair_geno == [3,2]:
        return 11

concor_table = defaultdict(lambda: [0 for x in range(0,12)])

while True:
    try:
        if query_pos == truth_pos:
            query_geno = geno_dict(query_line,query_id)
            truth_geno = geno_dict(truth_line,truth_id)
            #print(truth_geno)
            # iterate on the test set
            for key,value in query_geno.items():
                #print(value,truth_geno[key])
                # not considering missing in truth
                if not truth_geno[key] == 3:
                    match_index = concor_match(value,truth_geno[key])
                    concor_table[key][match_index] += 1
            query_line=next(query).strip().split()
            truth_line=next(truth).strip().split()
            query_pos=query_line[1]
            truth_pos=truth_line[1]
        elif query_pos < truth_pos:
            query_line=next(query).strip().split()
            query_pos=query_line[1]
        elif query_pos > truth_pos:
            truth_line=next(truth).strip().split()
            truth_pos=truth_line[1]
    # end of list, calculate stat
    except StopIteration:
        #print(concor_table)
        all_geno_match = [0 for x in range(0,12)]
        # per anim stat
        for key,value in concor_table.items():
            print(key,*value)
            for ind,x in enumerate(value):
                all_geno_match[ind] += x
        # combine stat
        print("All",*all_geno_match)
        # concordance stat
        print("Concordance", sum(all_geno_match[x] for x in [0,4,8])/sum(all_geno_match[x] for x in range(0,9)))
        print("NR Sensitivity", sum(all_geno_match[x] for x in [4,5,7,8])/sum(all_geno_match[x] for x in [1,2,4,5,7,8,10,11]))
        print("NR Discrepancy", sum(all_geno_match[x] for x in [2,3,5,6,7])/sum(all_geno_match[x] for x in range(1,9)))
        print("Precision", sum(all_geno_match[x] for x in [4,8])/sum(all_geno_match[x] for x in range(3,9)))
        break

    # Index of match query and truth as in the paper

        # truth
           # 0 1  2 
         # 0 0 1  2 
         # 1 3 4  5
    # test 2 6 7  8   
         # 3 9 10 11
         
    


