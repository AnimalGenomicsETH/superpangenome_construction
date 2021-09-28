#!/usr/bin/env python
assemb_dir = "/cluster/work/pausch/danang/psd/scratch/real_data"
assemb_list, = glob_wildcards(assemb_dir + "/assembly/{breed}_aut.fa")
# print(assemb_list)
ref_genome = "UCD"
chromo_list = list(range(1, 30))
# impgr_list=['hg','pg','og','vg','gfa','xg']
index_list=['xg','giraffe.gbwt','gg','min','snarl','dist']


#graph prog
# prog_list=["minigraph","pggb"]
# prog_list=["pggb",'minigraph','linear']
# prog_list=["pggb",'minigraph','cactus']
prog_list=["pggb",'minigraph','cactus']


# for pggb
dirwork = "/cluster/work/pausch/danang/psd/scratch/real_data"
sifdir = "/cluster/work/pausch/danang/psd/bin/sif"

#for mapping
anim_list=["BSWCHEM110294048847"]
fastq_dir="/cluster/work/pausch/inputs/fastq/BTA"


#toggle for optional analysis, change 0 to 1 to run respective analysis
run_mapping=1
run_cactus=1

def optional_output(run_mapping,run_cactus):
    opt_out=[]
    if run_mapping:
        #opt_out.extend(expand("mapped/{prog}_{anim}.gam",prog=prog_list,anim=anim_list)),
        opt_out.extend(expand("mapped/{prog}_{anim}_mapping_stat_up.tsv",prog=prog_list,anim=anim_list)),
        opt_out.extend(expand("mapped/{prog}_{anim}_up.gam",prog=prog_list,anim=anim_list))
    if run_cactus:
        opt_out.append("graph/cactus/cactus_combined.vg"),
        opt_out.append("graph/cactus/cactus_combined.gfa"),
        opt_out.append("graph/cactus/cactus_combined_chop.vg"),
        opt_out.append("graph/cactus/graph_cactus_chop.gfa"),
        #opt_out.append("combine_finished.tsv")
    return opt_out


rule all:
    input: expand("graph/minigraph/graph_{chromo}.gfa", chromo=chromo_list),
           expand("graph/pggb_{chromo}/pggb_{chromo}.gfa", chromo=chromo_list),
           #expand("graph/{prog}/graph_{prog}_test.{impgr}",prog=prog_list,impgr=impgr_list),
        #    expand("graph/{prog}/graph_{prog}_chop.vg",prog=prog_list)
           #expand("graph/pggb/graph_pggb_chop_test.{impgr}", impgr=impgr_list),
           expand("graph/{prog}/graph_{prog}.{ind}",prog=prog_list, ind=index_list),
           optional_output(run_mapping,run_cactus),
           expand("graph/{prog}/graph_{prog}_stat.tsv",prog=prog_list),
           

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
        full_graph_gfa="graph/minigraph/graph_minigraph_combined.gfa",
        chop_graph="graph/minigraph/graph_minigraph_chop.vg",
        chop_graph_gfa="graph/minigraph/graph_minigraph_chop.gfa"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """
        
        vg combine {input} > {output.full_graph}

        vg convert -f {output.full_graph} > {output.full_graph_gfa}

        vg mod -X 256 {output.full_graph} > {output.chop_graph}

        vg convert -t {threads} -f {output.chop_graph} > {output.chop_graph_gfa}
 

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
        # xg="graph/pggb_{chromo}/pggb_{chromo}.xg",
        # pg="graph/pggb_{chromo}/pggb_{chromo}.pg",
        # hg="graph/pggb_{chromo}/pggb_{chromo}.hg",
        # og="graph/pggb_{chromo}/pggb_{chromo}.og",
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    params:
    shell:
        """
        vg convert --threads {threads} -g {input} -v > {output.vg}

        # vg convert --threads {threads} -g {input} -x > {output.xg}

        # vg convert --threads {threads} -g {input} -p > {output.pg}

        # vg convert --threads {threads} -g {input} -a > {output.hg}

        # vg convert --threads {threads} -g {input} -o > {output.og}

        """


# rule combine_pggb:
#     input:  expand("graph/pggb_{chromo}/pggb_{chromo}.{{impgr}}",chromo=chromo_list,impgr=impgr_list)
#     output: 
#         full_graph="graph/pggb/graph_pggb_test.{impgr}"
#     threads: 10
#     resources:
#         mem_mb = 5000,
#         walltime = "04:00"
#     shell:
#         """
        
#         vg combine {input} > {output.full_graph}
       
#         """


# rule convert_pggb_gfa:
#         input: "graph/pggb/graph_pggb_test.gfa",
#         output: "graph/pggb/graph_pggb_combined.gfa"
#         threads: 1
#         resources:
#            mem_mb=1000 ,
#            walltime= "01:00"
#         shell:
#            """
#             ln -s {input} > {output}
#            """



rule combine_pggb:
    input:  expand("graph/pggb_{chromo}/pggb_{chromo}.vg",chromo=chromo_list)
    output: 
        full_graph="graph/pggb/graph_pggb_combined.vg"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
        
        vg combine {input} > {output.full_graph}
       
        """


rule convert_pggb_gfa:
        input: "graph/pggb/graph_pggb_combined.vg",
        output: "graph/pggb/graph_pggb_combined.gfa"
        threads: 10
        resources:
           mem_mb=2000 ,
           walltime= "01:00"
        shell:
           """

            vg convert -t {threads} -f {input} > {output} 

           """

rule chop_pggb:
    input:  "graph/pggb/graph_pggb_combined.vg"
    output: 
        chop_graph="graph/pggb/graph_pggb_chop.vg",
        gfa="graph/pggb/graph_pggb_chop.gfa"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
        
        vg mod -t {threads} -X 1000 {input} > {output.chop_graph}

        vg convert -t {threads} -f {input} > {output.gfa} 
       
        """


#cactus pipeline
include: "cactus_pipe.smk"       



## statistics of the graph 
rule calculate_graph_statistics:
        input:"graph/{prog}/graph_{prog}_combined.gfa"
        output:"graph/{prog}/graph_{prog}_stat.tsv"
        threads:10
        resources:
           mem_mb=2000 ,
           walltime= "04:00"
        params:
            grtype="{prog}"
        shell:
           """
            
            ./graph_stat.py  -g {input}  -o {output} -r {{1..29}}_UCD -t {params.grtype} 

           """

#indexing require for giraffe mapping 

#- xg : `vg convert -x -g gfa_fixed.gfa > testreg.xg`
#- gbwt (.giraffe.gbwt): `vg gbwt  -G gfa_fixed.gfa  --path-regex '(.*)_(.*)' --path-fields _CS -o testreg.gbwt`
#- gbwt graph (.gg): `vg gbwt -g testreg.gg -x testreg.xg testreg.gbwt` 
#- minimizer (.min): `vg minimizer -t 32 -g testreg.giraffe.gbwt -i testreg.min testreg.xg`
#- snarl (.snarl): `vg snarls -t 32 -T testreg.xg > testreg.snarls`
#- distance (.dist) : build from minimizer and snarl `vg index -t 32 testreg.xg -s testreg.snarls -j testreg.dist`

rule create_index_xg:
        input: "graph/{prog}/graph_{prog}_chop.gfa"
        output: "graph/{prog}/graph_{prog}.xg"
        threads: 32
        resources:
           mem_mb= 2000 ,
           walltime= "01:00"
        params:
        shell:
           """

           vg convert -t {threads} -x -g {input} > {output}

           """

rule create_index_gbwt:
        input: "graph/{prog}/graph_{prog}_chop.gfa"
        output: "graph/{prog}/graph_{prog}.giraffe.gbwt"
        threads: 32
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        shell:
           """

        #    vg gbwt -G {input} \
        #    --path-regex '(.*)_(.*)' --path-fields _CS \
        #    -o {output}

           vg gbwt --index-paths -x {input} -o {output}

           """

rule create_gbwt_graph:
        input: 
            xg="graph/{prog}/graph_{prog}.xg",
            gbwt="graph/{prog}/graph_{prog}.giraffe.gbwt"
        output: "graph/{prog}/graph_{prog}.gg"
        threads: 32
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        shell:
           """

           vg gbwt -g {output} -x {input.xg} {input.gbwt}

           """

rule create_index_minimizer:
        input:
            xg="graph/{prog}/graph_{prog}.xg",
            gbwt="graph/{prog}/graph_{prog}.giraffe.gbwt"
        output: "graph/{prog}/graph_{prog}.min"
        threads: 32
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        shell:
           """

           vg minimizer -t {threads} -g {input.gbwt} -i {output} {input.xg}

           """

rule identify_snarl:
        input: "graph/{prog}/graph_{prog}.xg"
        output: "graph/{prog}/graph_{prog}.snarl"
        threads: 32
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        shell:
           """

           vg snarls -t {threads} -T {input} > {output}

           """

rule create_index_distance:
        input: xg="graph/{prog}/graph_{prog}.xg",
               snarl="graph/{prog}/graph_{prog}.snarl"
        output:"graph/{prog}/graph_{prog}.dist"
        threads: 32
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        shell:
           """

           vg index -t {threads} {input.xg} \
            -s {input.snarl} -j {output}

           """

#mapping with giraffe 
rule giraffe_mapping:
        input:
            xg="graph/{prog}/graph_{prog}.xg",
            gg="graph/{prog}/graph_{prog}.gg",
            gbwt="graph/{prog}/graph_{prog}.giraffe.gbwt",
            mini="graph/{prog}/graph_{prog}.min",
            dist="graph/{prog}/graph_{prog}.dist",
            r1=fastq_dir + "/{anim}_R1.fastq.gz",
            r2=fastq_dir + "/{anim}_R2.fastq.gz"
        output:"mapped/{prog}_{anim}.gam"
        threads:32
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        params:
            fastq_dir=fastq_dir
        shell:
           """

           vg giraffe -t {threads} \
           -x {input.xg} \
           -g {input.gg} \
           -H {input.gbwt} \
           -m {input.mini} \
           -d {input.dist} \
           -p -f {input.r1} -f {input.r2} > {output}

           """


#updated indexing in version 1.34 
#this combine gbwt and gg into gbz file 

rule create_gbz_index:
        input: 
            xg="graph/{prog}/graph_{prog}.xg",
            gbwt="graph/{prog}/graph_{prog}.giraffe.gbwt"
        output:"graph/{prog}/graph_{prog}.giraffe.gbz"
        threads: 10
        resources:
           mem_mb= 5000 ,
           walltime= "04:00"
        shell:
           """

            vg gbwt --gbz-format -g {output} -x {input.xg} {input.gbwt}

           """

#mapping with giraffe 
rule giraffe_mapping_updated:
        input:
            xg="graph/{prog}/graph_{prog}.xg",
            gbz="graph/{prog}/graph_{prog}.giraffe.gbz",
            mini="graph/{prog}/graph_{prog}.min",
            dist="graph/{prog}/graph_{prog}.dist",
            r1=fastq_dir + "/{anim}_R1.fastq.gz",
            r2=fastq_dir + "/{anim}_R2.fastq.gz"
        output:"mapped/{prog}_{anim}_up.gam"
        threads:32
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        params:
            fastq_dir=fastq_dir
        shell:
           """

           vg giraffe -t {threads} \
           -x {input.xg} \
           -Z {input.gbz} \
           -m {input.mini} \
           -d {input.dist} \
           -p -f {input.r1} -f {input.r2} > {output}

           """

rule graph_mapping_statistics:
        input:"mapped/{prog}_{anim}_up.gam"
        output:"mapped/{prog}_{anim}_mapping_stat_up.tsv"
        threads:10
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        shell:
           """
        
            vg stats -a {input} > {output}
           
           """   


### conventional vg mapping 

rule prune_vg_graph:
        input:"graph/{prog}/graph_{prog}_chop.vg"
        output:"graph/{prog}/graph_{prog}_pruned.vg"
        threads:10
        resources:
           mem_mb=5000 ,
           walltime= "04:00"
        shell:
           """
            
            vg prune -t {threads} -M 32 --restore-paths {input} > {output} 

           """

rule index_vg_gcsa:
        input:"graph/{prog}/graph_{prog}_pruned.vg"
        output:"graph/{prog}/graph_{prog}.gcsa"
        threads:10
        resources:
           mem_mb=2000 ,
           walltime= "04:00"
        params:
           scrdir="graph/{prog}/graph_{prog}_pruned.vg"
        shell:
           """

           vg index --temp-dir {params.tmp_dir} -p -t {threads} -g {output.gcsa} {input} 

           """

rule vg_map_conventional:
        input:
            gcsa="graph/{prog}/graph_{prog}.gcsa",
            xg="graph/{prog}/graph_{prog}.xg",
            gbwt="graph/{prog}/graph_{prog}.giraffe.gbwt",
            r1=fastq_dir + "/{anim}_R1.fastq.gz",
            r2=fastq_dir + "/{anim}_R2.fastq.gz"
        output:"mapped/{prog}_{anim}_conv.gam"
        threads:32
        resources:
           mem_mb=2000 ,
           walltime= "04:00"
        shell:
           """
            
            vg map -t {threads} \
            -x {input.xg} -g {input.gcsa} -1 {index.gbwt} \
            -f {input.r1} -f {input.r2} > {output}

           """

rule graph_mapping_statistics_conv:
        input:"mapped/{prog}_{anim}_conv.gam"
        output:"mapped/{prog}_{anim}_mapping_stat_conv.tsv"
        threads:10
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        shell:
           """
        
            vg stats -a {input} > {output}
           
           """   
