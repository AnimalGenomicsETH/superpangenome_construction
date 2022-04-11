#!/usr/bin/env python

breeds_list=["BSW"]
fasta_dir="/cluster/work/pausch/danang/psd/scratch/real_comp/assembly"
snarl_dir="/cluster/work/pausch/danang/psd/scratch/real_comp/graph/snarl"
chromo_list=range(25,30)
prog_list=["cactus","pggb"]

rule all:
    input:
        expand("{breeds}/{breeds}_overlap_stat",breeds=breeds_list)

rule map_assembly:
    input:
        ref_genome= fasta_dir + "/{chromo}/UCD_{chromo}.fa",
        alt_genome = fasta_dir + "/{chromo}/{breeds}_{chromo}.fa"
    output: "{breeds}/{breeds}_{chromo}.paf"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "02:00"
    shell:
        """

        minimap2 -cx asm5 -t{threads} \
                --cs {input.ref_genome} {input.alt_genome} |
                sort -k6,6 -k8,8n > {output}
        """


rule call_variations:
    input:"{breeds}/{breeds}_{chromo}.paf"
    output:"{breeds}/{breeds}_{chromo}_varassemb.tsv"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "00:10"
    shell:
        """

        paftools.js call {input} > {output}

        """

rule overlap_variations:
    input:
        snarl_var = snarl_dir + "/{prog}/{prog}_{chromo}_snarl.vcf",
        paf_var = "{breeds}/{breeds}_{chromo}_varassemb.tsv"
    output: "{breeds}/{breeds}_{chromo}_{prog}_paf_overlap.tsv"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "00:30"
    shell:
        """
        
         ./con_dec.py {input.snarl_var} {wildcards.breeds} {input.paf_var} |
         awk -v prog={wildcards.prog} -v chromo={wildcards.chromo} '{{print prog,chromo,$0}}' OFS="\\t" > {output}

        """

rule combine_stat:
    input:expand("{{breeds}}/{{breeds}}_{chromo}_{prog}_paf_overlap.tsv",chromo=chromo_list,prog=prog_list)
    output:"{breeds}/{breeds}_overlap_stat"
    threads: 10
    resources:
        mem_mb= 1000 ,
        walltime= "00:10"
    shell:
        """
        
        cat {input} > {output}

        """



