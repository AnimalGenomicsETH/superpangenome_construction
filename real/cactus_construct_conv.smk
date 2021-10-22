#!/usr/bin/env python

dirwork="/cluster/work/pausch/danang/psd/scratch/real_yak2"
#chromo_list = list(range(25, 26))
chromo_list = [25]
#assemb_list = ["UCD","Angus","Highland","OBV","Brahman","Yak","BSW","Pied","Gaur","Nellore"]
assemb_list = ["UCD","Angus"]
ref_genome = "UCD"
rep_library="BosTau9_repeat_library.fasta"

rule all:
    input:
        expand("assembly/{chromo}/{chromo}_rep/{breed}_{chromo}.fa.masked",chromo=chromo_list,breed=assemb_list),
        expand("graph/cactus/cactus_{chromo}_seqfile.tsv",chromo=chromo_list)


    
rule cactus_masking_chromosome:
    input:"assembly/{chromo}/{breed}_{chromo}.fa"
    output:"assembly/{chromo}/{chromo}_rep/{breed}_{chromo}.fa.masked"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    params:
        rep_library=rep_library
    shell:
        """

        RepeatMasker -dir assembly/{wildcards.chromo}/{wildcards.chromo}_rep -no_is \
                -qq -xsmall -lib {params.rep_library} {input}

        """

rule sketch_assembly:
    input: "assembly/{breed}_aut.fa"
    output: "tree/{breed}.fa.msh"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        mash sketch -p {threads} -o {output} {input}

        """

rule combined_sketch:
    input: expand("tree/{breed}.fa.msh", breed=assemb_list)
    output: "tree/combined_sketch.msh"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        mash paste {output} {input}

        """

rule estimate_distance:
    input: "tree/combined_sketch.msh"
    output: "tree/combined_distance.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

        mash dist {input} {input} > {output}

        """


rule cactus_create_seqfile:
    input:
        assemb=expand("assembly/{{chromo}}/{{chromo}}_rep/{breed}_{{chromo}}.fa.masked",breed=assemb_list),
        distance_file=rules.estimate_distance.output
    output:"graph/cactus/cactus_{chromo}_seqfile.tsv"
    threads: 10
    resources:
        mem_mb=1000,
        walltime="01:00"
    params:
        fasta_dir="assembly/{chromo}/{chromo}_rep"
    shell:
        """

        ./construct_tree.R {input.distance_file} {params.fasta_dir} {output} fa.masked

        """


# rule cactus_align:


# rule cactus_combine_graph:



