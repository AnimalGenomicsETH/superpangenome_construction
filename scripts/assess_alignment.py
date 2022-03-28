#!/usr/bin/env python
"""

Assess alignment assemblies back into graphs
There are many hardcoded filepaths

"""

chromo_list=range(20,25)
chromo_dir="/cluster/work/pausch/danang/psd/scratch/real_comp/comb_chromo"


#get breed list
breed_list = [ x.strip() for x in open("breed_list.tsv").readlines() ]
program_list = ["cactus","pggb","minigraph"]

rule all:
    input:
        expand("stat_graph/{chromo}_{prog}_{anim}_pathstat.tsv",chromo=chromo_list,prog=program_list,anim=breed_list),
        expand("stat_graph/{chromo}_{prog}_{anim}_dist.tsv",chromo=chromo_list,prog=program_list,anim=breed_list)


rule split_fasta:
    input:chromo_dir + "/{chromo}_comb.fa"
    output:
        breed_fasta="splitfa/{chromo}/{chromo}_{anim}.fa",
        split_fasta="splitfa/{chromo}/split_{chromo}_{anim}.fa",
    threads:2
    resources:
        mem_mb=2000 ,
        walltime= "01:00"
    shell:
        """

        samtools faidx {input} {wildcards.chromo}_{wildcards.anim} > {output.breed_fasta}
        
        split_fa.py {output.breed_fasta} 1000000 > {output.split_fasta}

        """

def get_graph_file(graph,chromo):
    basedir="/cluster/work/pausch/danang/psd/scratch/real_comp/graph"
    if graph=="cactus":
        return basedir + f"/cactus/cactus_simple_{chromo}.gfa"
    elif graph=="pggb":
        return basedir + f"/pggb_{chromo}/pggb_{chromo}.gfa"   
    elif graph=="minigraph":
        return basedir + f"/minigraph/graph_{chromo}_path.gfa"


rule align_to_graph:
    input:
        fasta="splitfa/{chromo}/split_{chromo}_{anim}.fa",
        graph=lambda wildcards: get_graph_file(wildcards.prog,wildcards.chromo)
    output:"align_graph/{chromo}_{prog}_{anim}.gaf"
    threads:15
    resources:
        mem_mb= 5000,
        walltime= "01:00"
    shell:
        """

        module load gcc/6.3.0

        GraphAligner -g {input.graph} -f {input.fasta} --multimap-score-fraction 1 -a {output} -t 10 -x vg
        
        """

rule check_path:
    input:
        alignment="align_graph/{chromo}_{prog}_{anim}.gaf",
        graph=lambda wildcards: get_graph_file(wildcards.prog,wildcards.chromo)
    output:"stat_graph/{chromo}_{prog}_{anim}_pathstat.tsv"
    threads:10
    resources:
        mem_mb= 1000,
        walltime= "01:00"
    shell:
        """

        check_path.py -g {input.graph} -a {input.alignment} > {output}

        """

rule check_distance:
    input:"align_graph/{chromo}_{prog}_{anim}.gaf"
    output:"stat_graph/{chromo}_{prog}_{anim}_dist.tsv"
    threads:2
    resources:
        mem_mb=1000 ,
        walltime= "01:00"
    shell:
        """

        all_align=$(awk '{{cs+=$11}}END{{print cs}}' {input})
        edit_dis=$(awk '{{print $13}}' {input} | awk -F":" '{{cs+=$3}}END{{print cs}}')
        echo $all_align $edit_dis  $( echo "scale=2; $edit_dis*100/$all_align" | bc -l)

        """




