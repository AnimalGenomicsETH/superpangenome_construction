#!/usr/bin/env python


import subprocess 

breed_list="UCD,Angus,Highland,OBV,Brahman,Yak,BSW,Pied,Gaur,Nellore,Simmental,Bison"
breed_list=breed_list.split(",")

outgraph_dir="/cluster/scratch/cdanang/snarl_decomp"

rule all:
    input:
        expand("vnsplit/split_vntr_{chromo}.tsv",chromo=range(1,30)),
        expand(outgraph_dir + "/pggb_graph/pggb_{chromo}.og",chromo=range(1,30)),
        expand("vnsplit/vntr_coord_{chromo}_{breed}_pggb.tsv",chromo=range(1,3),breed=breed_list),
        expand("vnsplit/vntr_comb_{chromo}.tsv",chromo=range(1,3)),
        expand("vnsplit/{chromo}_vnseq_count.tsv",chromo=range(1,3))

rule split_vntr:
    input:"VNTRs.70.csv"
    output:"vnsplit/split_vntr_{chromo}.tsv"
    threads: 2
    resources:
        mem_mb= 2000,
        walltime= "02:00"
    shell:
        """

        awk -F"," -v chromo={wildcards.chromo} '$1 == chromo {{ print chromo"_UCD",$2,$3  }}' OFS="\\t" {input} > {output}
        
        """


rule build_odgi:
    input:"/cluster/work/pausch/danang/psd/scratch/real_comp/graph/pggb_{chromo}/pggb_{chromo}.gfa"
    output:outgraph_dir + "/pggb_graph/pggb_{chromo}.og"
    threads: 10
    resources: 
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        odgi build -g {input} -o {output} -t {threads}

        """

rule liftover_coordinate:
    input:
        chromo=outgraph_dir + "/pggb_graph/pggb_{chromo}.og",
        vntr="vnsplit/split_vntr_{chromo}.tsv"
    output:"vnsplit/vntr_coord_{chromo}_{breed}_pggb.tsv"
    threads: 20
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        odgi position -i {input.chromo} -b {input.vntr} -r {wildcards.chromo}_{wildcards.breed} -I > {output}


        """


rule concat_info:
    input:
        expand("vnsplit/vntr_coord_{{chromo}}_{breed}_pggb.tsv",breed=breed_list)
    output:
        "vnsplit/vntr_comb_{chromo}.tsv"
    threads: 5
    resources:
        mem_mb= 1000,
        walltime= "04:00"
    run:
        from collections import defaultdict
        coord_comb = defaultdict(dict)

        for input_file in input:
            with open(input_file) as infile:
                for line in infile:
                    token=line.strip().split()
                    chromo=token[0].split("_")[0]
                    breed=token[3].split(",")[0].split("_")[1]
                    vnid=f"{chromo}_{token[1]}_{token[2]}"
                    coord_comb[vnid][breed] = [chromo,breed,int(token[3].split(",")[1]),int(token[4].split(",")[1]),token[-1]]

        #order 
        outfile=open(output[0],"w")
        print("vnid",*breed_list,file=outfile)
        for vnid, value in coord_comb.items():
            print(vnid, *[",".join(str(x) for x in value.get(br,[br,0])) for br in breed_list],file=outfile)

        outfile.close()
        

fastadir="/cluster/work/pausch/danang/psd/scratch/assembly"
offset=500

def extract_fasta_region(line):
    sequence={}
    fastadir="/cluster/work/pausch/danang/psd/scratch/assembly"
    for comp in line:
        token=comp.split(",")
        if len(token) > 2:
            chromo,breed,start,stop,strand=comp.split(",")
            start, stop=int(start), int(stop)
            if strand == "+":
                seqfa=subprocess.run(f"samtools faidx {fastadir}/{chromo}/{breed}_{chromo}.fa {chromo}_{breed}:{start-offset}-{stop+offset} | seqtk2 seq -l 0",
                        shell=True,capture_output=True).stdout.decode("utf-8").split("\n")[1].upper()
                sequence[breed] = seqfa
            elif strand == "-":
                seqfa=subprocess.run(f"samtools faidx {fastadir}/{chromo}/{breed}_{chromo}.fa {chromo}_{breed}:{stop-offset}-{start+offset} | seqtk2 seq -l 0 -r",
                                    shell=True,capture_output=True).stdout.decode("utf-8").split("\n")[1].upper()
                sequence[breed] = seqfa
        else:
            sequence[token[0]] = "XX"
    return sequence
            

import regex
from collections import defaultdict

def count_VNTRs(sequences,TR):
    allowed_errors = int(len(TR)*0.2)
    counts = defaultdict(list)
    for asm, sequence in sequences.items():
        for hit in regex.finditer(f"(?eV1)({TR}){{e<={allowed_errors}}}",sequence,concurrent=True):
            counts[asm].append(sum(hit.fuzzy_counts))
    
    return counts

import math
from multiprocessing import Pool

trcomb = dict()
with open("VNTRs.70.csv") as infile:
    next(infile)
    for line in infile:
        comp=line.strip().split(",")
        trid=f"{comp[0]}_{comp[1]}_{comp[2]}"
        trcomb[trid] = comp[3]

from numpy import var as variance

def processs_line(line):
    token=line.strip().split()
    vnid=token[0]
    trseq=trcomb[vnid]
    combseq=extract_fasta_region(token[1:])

    counts = count_VNTRs(combseq,trseq)
    # return counts 

    stats = []
    for B in breed_list:
        if B in counts:
            stats.extend([len(counts[B]),sum(counts[B]),f'{variance(counts[B]):.3f}'])
        else:
            stats.extend([math.nan]*3)

    return ','.join(map(str, [vnid] + stats))
    # return 

from itertools import product

rule detect_pattern:
    input:
        vntrfile="VNTRs.70.csv",
        regfile="vnsplit/vntr_comb_{chromo}.tsv",
        fastafile= expand(fastadir + "/{{chromo}}/{breed}_{{chromo}}.fa",breed=breed_list)
    output:
        "vnsplit/{chromo}_vnseq_count.tsv"
    threads: 20
    resources:
        mem_mb= 1000,
        walltime= "04:00"
    run:
        with open(input.regfile) as infile, open(output[0],"w") as outfile:
            anims=next(infile).strip().split()[1:]
            print("vnid," + ','.join('_'.join(P) for P in product(breed_list,['count','sum','var'])),file=outfile)
            with Pool(18) as p:
                for result in p.imap_unordered(processs_line,infile,18):
                    print(result,file=outfile)
