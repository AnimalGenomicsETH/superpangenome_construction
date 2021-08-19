#!/usr/bin/env python

assemb_dir = "/cluster/work/pausch/danang/psd/scratch/real_data"
assemb_list, = glob_wildcards(assemb_dir + "/assembly/{breed}_aut.fa")

ref_genome = "UCD"
chromo_list = list(range(1, 30))

dirwork = "/cluster/work/pausch/danang/psd/scratch/sim/sim_rep2"
scrdir="/cluster/scratch/cdanang/scr_cac"

#make it chromosome-by-chromosome 
#this use input GFA from minigraph, aligned, and masked in a single pre-process step 


rule all:
    input: expand("graph/cactus/{chromo}/seqfile_{chromo}.txt",chromo=chromo_list),
           expand("graph/cactus/{chromo}/seqfile_{chromo}_masked.txt",chromo=chromo_list),
           expand("graph/cactus/{chromo}/cactus_{chromo}.paf",chromo=chromo_list),
           expand("graph/cactus/{chromo}/{anim}_{chromo}_masked.fa",anim=assemb_list,chromo=chromo_list)

localrules: create_seq_file 
rule create_seq_file:
        input:
            graph="graph/minigraph/graph_{chromo}.gfa",
            fasta=expand("assembly/{{chromo}}/{anim}_{{chromo}}.fa",anim=assemb_list)
        output:
            seqfile="graph/cactus/{chromo}/seqfile_{chromo}.txt",
            seqfile_masked="graph/cactus/{chromo}/seqfile_{chromo}_masked.txt"
        params:
            prefix="graph/cactus/{chromo}",
            anims=assemb_list,
            chromo="{chromo}"
        shell:
           """

           for anim in {params.anims}
           do
                echo {params.chromo}_${{anim}} assembly/{params.chromo}/${{anim}}_{params.chromo}.fa >> {params.prefix}/seqfile_{params.chromo}.txt
                echo {params.chromo}_${{anim}} {params.prefix}/${{anim}}_{params.chromo}_masked.fa >> {params.prefix}/seqfile_{params.chromo}_masked.txt
           done 

           """

rule cactus_preprocess_map:
        input: 
            graph="graph/minigraph/graph_{chromo}.gfa",
            fasta=expand("assembly/{{chromo}}/{anim}_{{chromo}}.fa",anim=assemb_list),
            seqfile="graph/cactus/{chromo}/seqfile_{chromo}.txt",
        output: 
            paf="graph/cactus/{chromo}/cactus_{chromo}.paf"
        threads: 10
        resources:
           mem_mb= 5000,
           walltime= "04:00"
        params:
            prefix="graph/cactus/{chromo}",
            scrdir=scrdir +"/jobstore_map_{chromo}",
            anims=assemb_list,
            chromo="{chromo}"
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH


            cactus-graphmap $PWD/{params.chromo}_src  \
            {input.seqfile} \
            {input.graph} {output.paf} \
            --outputFasta {params.prefix}/cactus_{wildcards.chromo}.fa \
            --realTimeLogging

           """

rule cactus_preprocess_masked:
        input:
            seqfile="graph/cactus/{chromo}/seqfile_{chromo}.txt",
            seqfile_masked="graph/cactus/{chromo}/seqfile_{chromo}_masked.txt"
        output:expand("graph/cactus/{{chromo}}/{anim}_{{chromo}}_masked.fa",anim=assemb_list)
        threads:32 
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        params:
            chromo = "{chromo}"
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH


           cactus-preprocess $PWD/{params.chromo}_src_masked \
           {input.seqfile} {input.seqfile_masked} \
           --realTimeLogging --binariesMode local \
           --maskAlpha --brnnCores 8

           """

