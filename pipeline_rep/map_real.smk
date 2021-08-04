#!/usr/bin/env python 
assemb_dir="/cluster/work/pausch/danang/psd/scratch/real_data"
assemb_list,=glob_wildcards(assemb_dir+"/assembly/{breed}_aut.fa")
# print(assemb_list)
ref_genome="UCD"
chromo_list=list(range(25,27))

rule all:
    input : expand("graph/minigraph/graph_{chromo}.gfa",chromo=chromo_list)

rule sketch_assembly:
    input : "assembly/{assemb}_aut.fa"
    output: "tree/{assemb}.fa.msh"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        mash sketch -p {threads} -o {output} {input}

        """

localrules: combined_sketch
rule combined_sketch:
    input: expand("tree/{assemb}.fa.msh", assemb=assemb_list)
    output: "tree/combined_sketch.msh"
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

localrules: split_chromosome
rule split_chromosome:
    input:"assembly/{assemb}_aut.fa"
    output:"assembly/{chromo}/{assemb}_{chromo}.fa"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    params: ref = ref_genome
    shell:
        """

        if [[ "{wildcards.assemb}" == "{params.ref}" ]]; then 
            samtools faidx {input} {wildcards.chromo} > {output}
        else
            samtools faidx {input} {wildcards.chromo}_{wildcards.assemb} > {output}
        fi 

        """

localrules: determine_order
rule determine_order:
    input:"tree/combined_distance.tsv"
    output:"breed_order.tsv"
    params:
        ref=ref_genome
    run:
        import pandas as pd 

        allin=pd.read_csv(input[0],sep="\t",header=None)
        allin.columns=["anim1","anim2","dist","comp1","comp2"]
        allin["anim1"] = allin.anim1.str.extract(r"([A-Z][^_]+)")
        allin["anim2"] = allin.anim2.str.extract(r"([A-Z][^_]+)")

        allsel=allin.loc[allin.anim1==params.ref]


        allsort=allsel.sort_values(by="dist")
        outfile=open(output[0],"w")
        print(*list(allsort.anim2),file=outfile)

def get_order_list(order_file,chromo):
    infile=open(order_file)
    return [f"assembly/{chromo}/{br}_{chromo}.fa" for br in infile.readlines()[0].strip().split(" ")]
    infile.close()

rule run_minigraph:
    input:
        order_file="breed_order.tsv",
        assembly=expand("assembly/{{chromo}}/{breed}_{{chromo}}.fa",breed=assemb_list)
    output:"graph/minigraph/graph_{chromo}.gfa"
    threads: 10
    resources:
       mem_mb= 2000,
       walltime= "04:00"
    params: order_list=get_order_list("breed_order.tsv","{chromo}")
    shell:
        """

        # cat minigraph {params.order_list} > {output}

        minigraph -t {threads} -xggs {params.order_list} > {output}

        """