#!/usr/bin/env python

import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description = __doc__,
                    formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i","--input",help="input file"),
    parser.add_argument("-o","--output",help="output file")
    parser.add_argument("-l","--statlen",action="store_true",help="Use this tag to print snarl len (stdout)")
    parser.add_argument("-v","--svonly",action="store_true",help="Use this to extract sv snarls (stdout)")
    return parser.parse_args()


if __name__ == "__main__":
    args=parse_args()
    input_file=args.input
    statlen=args.statlen
    all_var=0
    max_len=0
    snp_var=0
    bial_var=0
    multi_var=0
    sv_var=0
    subt=0
    ins=0
    delt=0
    subt_sv=0
    ins_sv=0
    delt_sv=0
    with open(input_file) as infile:
        for line in infile:
            if line.startswith("#CHROM") and (statlen or args.svonly):
                if args.svonly:
                    line_comp=line.strip().split()
                    print(*line_comp[0:4], *[x.split(".")[0] for x in line_comp[9:len(line_comp)]],sep="\t")
                else:
                    print(line.strip()
            if not line.startswith('#'):
                all_var+=1
                token=line.strip().split()
                ref=token[3]
                ref_len=len(ref)
                alt=token[4].split(",")
                alt_len=[len(x) for x in alt]
                all_len=alt_len + [ref_len]
                if len(alt) > 1:
                    multi_var += 1
                if len(alt) == 1:
                    bial_var += 1
                if max(all_len) >= 50:
                    sv_var += 1
                    if args.svonly:
                        print(*token,sep="\t")
                        continue
                if max(all_len) < 50:
                    snp_var += 1
                if max(all_len) > max_len:
                    max_len = max(all_len)
                for comp in alt_len:
                    if comp == ref_len:
                        subt += 1
                    if comp > ref_len:
                        ins += 1
                    if comp < ref_len:
                        delt += 1
                    if max(all_len) >= 50:
                        if comp == ref_len:
                            subt_sv += 1
                        if comp > ref_len:
                            ins_sv += 1
                        if comp < ref_len:
                            delt_sv += 1
                if statlen:
                    alt_len_list=[]
                    for x in token[9:len(token)]:
                        if re.search(r"/|\.",x):
                            alt_len_list.append(x)
                        else:
                            x=int(x)
                            if x:
                                alt_len_list.append(alt_len[int(x)-1])
                            else:
                                alt_len_list.append('R')
                    #print(token[0],token[1],token[2],ref_len,alt_len,*token[9:len(token)],sep="\t")
                    print(token[0],token[1],token[2],ref_len,*alt_len_list,sep="\t")
    if not statlen and not args.svonly:
        print(all_var,max_len,snp_var,sv_var,bial_var,multi_var,subt,ins,delt,subt_sv,ins_sv,delt_sv)






