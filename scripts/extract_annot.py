#!/usr/bin/env python
"""

Take ensemble annot VCF  and simplify the mutation report

Input: Jasmine VCF with annotation
Output: Stdout with tab |chromo|location|merged sample|mutation loc (coding, intergenic, intronic) | mutation effects (modifier,low,medium,high)

"""

import sys
import re

with open(sys.argv[1]) as infile:

    effect_rank = {"MODIFIER":0,
                   "LOW":1,
                   "MODERATE":2,
                   "HIGH":3}

    inv_effect = {v:k for k,v in effect_rank.items()}

    for line in infile:
        if not line.startswith("#"):
            token=line.strip().split()
            suppvec = re.findall("SUPP_VEC=([^;]+);",token[7])[0]
            alt = token[4].split(",")
            annot = token[7].split("|")
            ind_vartype = [1]
            while True:
                if ind_vartype[-1]+26 < len(annot)-1:
                    ind_vartype.append(ind_vartype[-1]+26)
                else:
                    break
            #print(len(annot))
            #print(ind_vartype)
            all_effect = [annot[x+1] for x in ind_vartype]
            max_effect = 0
            for comp in all_effect:
                if effect_rank[comp] > max_effect:
                    max_effect = effect_rank[comp]
            max_effect = inv_effect[max_effect]
            symbol_id = [annot[x+2] for x in ind_vartype]
            symbol_id = set([x if x else "n" for x in symbol_id])
            gene_id = [annot[x+3] for x in ind_vartype]
            gene_id = set([x if x else "n" for x in gene_id])
            print(token[0], token[1], suppvec, 
                set([annot[x] for x in ind_vartype]),
                set([annot[x+1] for x in ind_vartype]), 
                ",".join(list(symbol_id)),
                ",".join(list(gene_id)),
                max_effect)
