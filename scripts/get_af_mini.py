#!/usr/bin/env python
from collections import Counter

#extract af from combine minigraph

with open("sv_call_node_mini.bed") as infile:
    infoline=next(infile).strip().split(" ")
    anim=infoline[5:]
    for line in infile:
        token=line.strip().split()
        allele_list=[x.split(":")[0] for x in token[5:]]    
        allele_id=list(set(allele_list))
        ref_allele=allele_list[0]
        alt_allele=[x for x in allele_id if x != "." and x != ref_allele]
        alt_ind={x:ind+1 for ind,x in enumerate(alt_allele)}
        #decode the allele 
        decode_allele=[]
        for x in allele_list:
            if x == ".":
                decode_allele.append("NA")
            else:
                if x == ref_allele:
                    decode_allele.append(0)
                else:
                    decode_allele.append(alt_ind[x])
        missing_rate=sum(x=="." for x in allele_list)
        total_geno=len(allele_list) - missing_rate
        multial=1 if len(alt_ind) > 1 else 0 
        allele_counter=Counter(decode_allele)
        
        af_ab=[]
        af_rel=[]
        all_geno=sum(allele_counter.values())
        typed_geno=sum(x for key,x in allele_counter.items() if key != "NA")
        for key,value in allele_counter.items():
            if key != "NA" and key != 0:
                af_ab.append(value/all_geno)
                af_rel.append(value/typed_geno)
                
        
        #print(*token[:5],*allele_list,ref_allele,alt_allele)
        print(*token[:5],*decode_allele,af_ab,af_rel,missing_rate)
        #print(allele_counter)


