#!/usr/bin/env python
assemb_dir = "/cluster/work/pausch/danang/psd/scratch/real_data"
assemb_list, = glob_wildcards(assemb_dir + "/assembly/{breed}_aut.fa")
# print(assemb_list)
ref_genome = "UCD"
chromo_list = list(range(1, 30))
impgr_list=['hg','pg','og','vg','gfa','xg']


#graph prog
# prog_list=["minigraph","pggb"]
prog_list=["pggb"]

# for pggb
dirwork = "/cluster/work/pausch/danang/psd/scratch/real_data"
sifdir = "/cluster/work/pausch/danang/psd/bin/sif"

rule all:
    input: expand("graph/minigraph/graph_{chromo}.gfa", chromo=chromo_list),
           expand("graph/pggb_{chromo}/pggb_{chromo}.gfa", chromo=chromo_list),
           expand("graph/{prog}/graph_{prog}_test.{impgr}",prog=prog_list,impgr=impgr_list),
        #    expand("graph/{prog}/graph_{prog}_chop.vg",prog=prog_list)
           expand("graph/pggb/graph_pggb_chop_test.{impgr}", impgr=impgr_list),
            

rule sketch_assembly:
    input: "assembly/{assemb}_aut.fa"
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
    input: "assembly/{assemb}_aut.fa"
    output: "assembly/{chromo}/{assemb}_{chromo}.fa"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
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
    input: "tree/combined_distance.tsv"
    output: "breed_order.tsv"
    params:
        ref = ref_genome
    run:
        import pandas as pd

        allin = pd.read_csv(input[0], sep="\t", header=None)
        allin.columns = ["anim1", "anim2", "dist", "comp1", "comp2"]
        allin["anim1"] = allin.anim1.str.extract(r"([A-Z][^_]+)")
        allin["anim2"] = allin.anim2.str.extract(r"([A-Z][^_]+)")

        allsel = allin.loc[allin.anim1 == params.ref]

        allsort = allsel.sort_values(by="dist")
        outfile = open(output[0], "w")
        print(*list(allsort.anim2), file=outfile)


def get_order_list(order_file, chromo):
    infile = open(order_file)
    return [f"assembly/{chromo}/{br}_{chromo}.fa" for br in infile.readlines()[0].strip().split(" ")]
    infile.close()


rule run_minigraph:
    input:
        order_file = "breed_order.tsv",
        assembly = expand("assembly/{{chromo}}/{breed}_{{chromo}}.fa", breed=assemb_list)
    output: "graph/minigraph/graph_{chromo}.gfa"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    params: order_list = get_order_list("breed_order.tsv", "{chromo}")
    shell:
        """

        # cat minigraph {params.order_list} > {output}

        minigraph -t {threads} -xggs {params.order_list} > {output}

        """


localrules: modify_minigraph
rule modify_minigraph:
    input:  "graph/minigraph/graph_{chromo}.gfa"
    output: "graph/minigraph/graph_{chromo}_path.vg"
    shell:
        """
        vg convert -r 0 -g {input} -v > {output}

        """

rule combine_minigraph:
    input:  expand("graph/minigraph/graph_{chromo}_path.vg",chromo=chromo_list)
    output: 
        full_graph="graph/minigraph/graph_minigraph.vg",
        chop_graph="graph/minigraph/graph_minigraph_chop.vg"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """
        
        vg combine {input} > {output.full_graph}

        vg mod -X 256 {output.full_graph} > {output.chop_graph}

        """

localrules: combine_genome
rule combine_genome:
    input: expand("assembly/{{chromo}}/{breed}_{{chromo}}.fa", breed=assemb_list)
    output: "comb_chromo/{chromo}_comb.fa"
    shell:
        """

        cat {input} > {output}

        """

rule construct_pggb:
    input:
        "comb_chromo/{chromo}_comb.fa"
    output:
        "graph/pggb_{chromo}/pggb_{chromo}.gfa"
    threads: 32
    resources:
        mem_mb = 2000,
        walltime = "24:00"
    params:
        sifdir = sifdir,
        dirwork = dirwork
    shell:
        """

        singularity run --bind {params.dirwork} {params.sifdir}/pggb.sif \
        'pggb -i {input} -t {threads} -s 100000 -p 90 -n 10 \
        -S -o graph/pggb_{wildcards.chromo}'

        ln -s --relative graph/pggb_{wildcards.chromo}/*.smooth.gfa {output}

        """

rule convert_ppgb:
    input: "graph/pggb_{chromo}/pggb_{chromo}.gfa"
    output: 
        vg="graph/pggb_{chromo}/pggb_{chromo}.vg",
        xg="graph/pggb_{chromo}/pggb_{chromo}.xg",
        pg="graph/pggb_{chromo}/pggb_{chromo}.pg",
        hg="graph/pggb_{chromo}/pggb_{chromo}.hg",
        og="graph/pggb_{chromo}/pggb_{chromo}.og",
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    params:
    shell:
        """
        vg convert --threads {threads} -g {input} -v > {output.vg}

        vg convert --threads {threads} -g {input} -x > {output.xg}

        vg convert --threads {threads} -g {input} -p > {output.pg}

        vg convert --threads {threads} -g {input} -a > {output.hg}

        vg convert --threads {threads} -g {input} -o > {output.og}

        """

rule combine_pggb:
    input:  expand("graph/pggb_{chromo}/pggb_{chromo}.{{impgr}}",chromo=chromo_list,impgr=impgr_list)
    output: 
        full_graph="graph/pggb/graph_pggb_test.{impgr}"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
        
        vg combine {input} > {output.full_graph}
       
        """


rule chop_pggb:
    input:  "graph/pggb/graph_pggb_test.{impgr}"
    output: 
        chop_graph="graph/pggb/graph_pggb_chop_test.{impgr}"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
        
        vg mod -t {threads} -X 1000 {input} > {output.chop_graph}
       
        """