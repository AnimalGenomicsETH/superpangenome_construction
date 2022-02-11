#!/usr/bin/env python
"""

This will take two vcf and report the sites which are concordance between to stdout

The usage 

./concor_graph_assemb.py decons_vcf sample paf_bed

"""

import sys
import re
import gzip


def extract_decons_geno(input_file,sample_id="default",vcf_call=0):
    decon_list = dict()
    anim_index = 0
    #if input_file[-2:] == "gz":
    #    infile = gzip.open(input_file)
    infile = open(input_file)

    for line in infile:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                if vcf_call:
                    anim_index = 9
                else:
                    for index,anim in enumerate(line.strip().split()):
                        if re.search(sample_id,anim):
                            anim_index = index
        else:
            line_comp = line.strip().split()
            pos = line_comp[1]
            ref = line_comp[3]
            alt = {x+1:y for x,y in enumerate(line_comp[4].split(","))}
            geno = line_comp[anim_index].split(":")[0]
            if not re.search("0",geno) and not re.search("\.",geno):
                    geno_list = list(set(geno.split("/")))
                    if len(geno_list) == 1:
                        allele = alt[int(geno_list[0])]
                        decon_list[int(pos)] = [ref,allele]
    infile.close()
    return decon_list

def extract_paf_geno(input_file):
    paf_list = dict()
    with open(input_file) as infile:
        for line in infile:
            if line.startswith("V"):
                token=line.strip().split()
                pos1 = token[2]
                pos2 = token[3]
                ref = token[6]
                alt = token[7]
                if "-" not in ref:
                    ref = ref.upper()
                if "-" not in alt:
                    alt = alt.upper()
                paf_list[int(pos1)] = [ref,alt]
                paf_list[int(pos2)] = [ref,alt]
    return paf_list

if __name__ == "__main__":
    decons_geno = extract_decons_geno(sys.argv[1],sample_id=sys.argv[2])
    call_geno = extract_decons_geno(sys.argv[3],vcf_call=1)
    
    for key,value in decons_geno.items():
        if call_geno.get(key,0):
            print(key,*value,key,*call_geno[key])
    
            

