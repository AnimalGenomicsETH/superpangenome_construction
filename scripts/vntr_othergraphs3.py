#!/usr/bin/env python

import pandas as pd
import sys
import regex
from multiprocessing import Pool 
from collections import defaultdict


# intersection of pggb and cactus
# cat <(awk 'NR==1{print "chromo","startbub","endbub","nodeid",$0}' OFS="\t" VNTR_inter.bed ) \
# <(bedtools intersect -wo -a pggb_snarl_combloc.bed -b VNTR_inter.bed) > pggb_vntr_inter.bed 

# cat <(awk 'NR==1{print "chromo","startbub","endbub","nodeid",$0}' OFS="\t" VNTR_inter.bed ) \
# <(bedtools intersect -wo -a cactus_snarl_combloc.bed -b VNTR_inter.bed) > cactus_vntr_inter.bed

# full vntr info 

# datvn = pd.read_csv("VNTRs.912.csv")
# datvn["vnid"] = ["VN" + str(x) for x in range(1,datvn.shape[0]+1)]

# intersection of the vntr between graph 


def get_inter(interfile):
    datinter=pd.read_csv(interfile,sep="\t",header=0)

    # select the most representative overlap
    chromo_id=int(regex.findall(r"\d+",sys.argv[1])[0])
    datsel = datinter.sort_values("interlen",ascending=False).groupby("vnid").first().reset_index()
    datsel = datsel.loc[datsel.chromo==chromo_id,]

    return datsel 


def get_seq(vargraph,datsel):
    vntr_graph_seq=defaultdict(dict)
    with open(vargraph) as infile:
        for line in infile:
            if line[:6] == "#CHROM":
                anims=line.strip().split()[4:]
            if line[0] != "#":
                token=line.strip().split()
                allele = []
                if token[2] in datsel.nodeid.values:
                    allele.append(token[3])
                    allele.extend(token[4].split(","))
                    vnid=datsel.loc[datsel.nodeid==token[2],"vnid"].values[0]
                    vnseq=datsel.loc[datsel.vnid==vnid,"TR"].values[0]
                    #combseq = {}
                    vntr_graph_seq[vnid]["UCD"] = allele[0]
                    for anim,alid in zip(anims,token[9:]):
                        if "/" in alid:
                            vntr_graph_seq[vnid][anim] = allele[int(alid.split("/")[0])]
                        elif alid == ".":
                            vntr_graph_seq[vnid][anim] ="XX"
                        else:
                            vntr_graph_seq[vnid][anim] = allele[int(alid)]
    return vntr_graph_seq,anims

def count_VNTRs(sequences,TR):
    allowed_errors = int(len(TR)*0.2)
    return {asm: len(regex.findall(f"({TR}){{e<={allowed_errors}}}",sequence,concurrent=True)) for asm, sequence in sequences.items()}

def count_pattern(line):
    vnid=line[0]
    vnseq=datsel.loc[datsel.vnid==vnid,"TR"].values[0]
    combseq = vntr_graph_seq[vnid]
    if combseq:
        return [vnid,count_VNTRs(combseq,vnseq)]
    return [0,0]


if __name__ == "__main__":
    
    datsel=get_inter("cactus_vntr_inter.bed")
    chromo_basepath="/cluster/work/pausch/danang/psd/scratch/real_comp/graph/snarl/cactus/"
    vargraph = chromo_basepath + sys.argv[1]
    vntr_graph_seq,anims=get_seq(vargraph,datsel)
    anims.append("UCD")
    chromo_id=int(regex.findall(r"\d+",sys.argv[1])[0])
    
    # cactus 
    # with open(sys.argv[3], "w") as outfile:
        # if chromo_id == 1:
            # print("prog","vntr_id","chromo","startpos","endpos",*anims,"vntr_seq",sep="\t",file=outfile)
        # with Pool(10) as p:
            # for result in p.imap_unordered(count_pattern, datsel.values, 10):
                # if result[1]:
                    # vn_id = result[0]
                    # minirep=datsel.loc[datsel.vnid==vn_id,anims].values
                    # minirep=[int(x) for x in minirep[0]]
                    # cacrep=[result[1][anim] for anim in anims] 
                    # vnseq=datsel.loc[datsel.vnid==vn_id,"TR"].values[0]
                    # vnchromo=datsel.loc[datsel.vnid==vn_id,"chr"].values[0]
                    # startpos=datsel.loc[datsel.vnid==vn_id,"start"].values[0]
                    # endpos=datsel.loc[datsel.vnid==vn_id,"end"].values[0]
                    # print("mini",vn_id,vnchromo,startpos,endpos,*minirep,vnseq,sep="\t",file=outfile)
                    # print("cactus",vn_id,vnchromo,startpos,endpos,*cacrep,vnseq,sep="\t",file=outfile)
    
    datsel=get_inter("pggb_vntr_inter.bed")
    chromo_basepath="/cluster/work/pausch/danang/psd/scratch/real_comp/graph/snarl/pggb/"
    vargraph = chromo_basepath + sys.argv[2]
    vntr_graph_seq,anims=get_seq(vargraph,datsel)
    anims.append("UCD")
    chromo_id=int(regex.findall(r"\d+",sys.argv[1])[0])

    # pggb
    with open(sys.argv[3], "a") as outfile:
        with Pool(10) as p:
            for result in p.imap_unordered(count_pattern, datsel.values, 10):
                if result[1]:
                    vn_id = result[0]
                    cacrep=[result[1][anim] for anim in anims] 
                    vnseq=datsel.loc[datsel.vnid==vn_id,"TR"].values[0]
                    vnchromo=datsel.loc[datsel.vnid==vn_id,"chr"].values[0]
                    startpos=datsel.loc[datsel.vnid==vn_id,"start"].values[0]
                    endpos=datsel.loc[datsel.vnid==vn_id,"end"].values[0]
                    print("pggb",vn_id,vnchromo,startpos,endpos,*cacrep,vnseq,sep="\t",file=outfile)
