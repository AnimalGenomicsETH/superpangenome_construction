#!/usr/bin/env python
"""
Intersect sv breakpoint with GFF annotations 
and report the most specific features from each breakpoints

intergenic --> intronic ---> exonic ---> CDS

"""
import argparse
import re
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description = __doc__,
                    formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--input",help="input bed intersect file")
    parser.add_argument("-o","--output",help="output breakpoint file")
    return parser.parse_args() 

def get_feature_table(mergedbreak):
    """
    parse the gff line to get name and ID from each feature
    """
    for feature in mergedbreak:
        annot = [x.split("=") for x in feature.fields[-1].split(";")]
        try:
            annot2 = {x1: x2 for x1, x2 in annot}
            ID = annot2.get("ID", "noid")
            if ID not in "noid":
                ID = ID.split(":")[1]
            Name = annot2.get("Name", "noname")
            # chromo, start_break, start_node, stop_node, svid, featsel
            yield [feature[0], feature[1], feature[3], feature[4], feature[5], feature[8], ID, Name]
        except:
            continue

def get_gene_id(lastgff):
    """
    parse the gff last line to get name and ID from each feature
    """
    annot = [x.split("=") for x in lastgff.split(";")]
    try:
        annot2 = {x1: x2 for x1, x2 in annot}
        # print(annot2)
        ID = annot2.get("ID", "noid")
        if ID not in "noid":
            ID = ID.split(":")[1]
        Name = annot2.get("Name", "noname")
        return(ID,Name)
    except:
        return("noid","noname2")

def extract_important_feature(prevfeat, curfeat):
    """

    Return the most important feature from gff
    """
    priority = {"CDS": 4, "exon": 3, "mRNA": 2, "gene": 1}
    prevprior = priority.get(prevfeat, 0)
    curprior = priority.get(curfeat, 0)

    if curprior == 0 and prevprior == 0:
        return "intergenic"
    if curprior >= prevprior:
        return curfeat
    else:
        return prevfeat


def get_feat_sel(svtype, featid, featname, featsel):
    """
        Extract gene name from gff features 
    """
    if re.search(r"gene", svtype):
        # report gene name which likely to be the shortest
        featsel = featid if len(featname) > len(featid) else featname
        # if no name, report the gene ID
        if featsel == "noname":
            featsel = featid
    return featsel

if __name__ == "__main__":

    args=parse_args()
    previd=""
    prevfeat = ""
    curfeat=""
    prev_gene_id =""
    prev_gene_name=""
    prev_start=""
    prev_chromo=""
    svid=""

    with open(args.input) as infile, open(args.output,"w") as outfile:
        for index,line in enumerate(infile):
            line_comp = line.strip().split("\t")
            chromo,start_pos,stop_pos,svid,ensmbl,chromo2,svtype,*_,svcomp=line_comp
            gene_id,gene_name =get_gene_id(svcomp)

            if svid == previd:
                curfeat = extract_important_feature(prevfeat,svtype)
                if svtype == "gene":
                    prev_gene_id = gene_id
                    prev_gene_name = gene_name
                prevfeat = curfeat
                previd = svid
                prev_start = start_pos
                prev_chromo = chromo
            else:
                print(previd,prev_chromo,prev_start,prevfeat,prev_gene_id,prev_gene_name, file=outfile)
                prevfeat = "intergenic" if svtype =="." else svtype
                previd = svid
                prev_start = start_pos
                prev_gene_id = gene_id
                prev_gene_name = gene_name
                prev_chromo = chromo



