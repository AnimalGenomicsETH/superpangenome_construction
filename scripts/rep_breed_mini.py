#!/usr/bin/env python


breed_list="UCD,Angus,Highland,OBV,Brahman,Yak,BSW,Pied,Gaur,Nellore,Simmental,Bison"
breed_list=breed_list.split(",")

chromo=range(1,30)
basedir="/cluster/work/pausch/danang/psd/scratch/real_comp/graph/minigraphbs"

rule all:
    input:
        expand("nonref_mini/{breed}_{chromo}.fa",breed=breed_list,chromo=chromo_list)

rule extract_nonref:
    input:
        graph=basedir + "/minigraph_{chromo}_bs.gfa",
        callsv=basedir + "/call_sv/sv_{chromo}_{breed}_bs.bed"
    output:
        "nonref_mini/{chromo}_{breed}_nonref.fa"
    threads:2
    resources:
        mem_mb= 4000,
        walltime= "01:00"
    shell:
        """

         get_perbreed_nonref.py {input.graph} {input.callsv} > {output}

        """

