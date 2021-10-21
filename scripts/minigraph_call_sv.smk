#!/usr/bin/env python
assemb_dir = "/cluster/work/pausch/danang/psd/scratch/real_all"
#assemb_list, = glob_wildcards(assemb_dir + "/assembly/{breed}_aut.fa")
assemb_list = ["UCD","Angus","Highland","OBV","Brahman","Yak","BSW","Pied","Gaur","Nellore"]
# print(assemb_list)
ref_genome = "UCD"
chromo_list = list(range(1, 30))
#gff annotation file
GFF_file="bos_taurus_annot4.gff3"

rule all:
        input:
            "graph/minigraph/call_sv/sv_all_call.bed",
            "graph/minigraph/call_sv/sv_annot_combined.bed"



rule call_sv_minigraph:
        input:
            graph="graph/minigraph/graph_{chromo}.gfa",
            assembly="assembly/{chromo}/{breed}_{chromo}.fa"
        output:"graph/minigraph/call_sv/sv_{chromo}_{breed}.bed"
        threads:10
        resources:
           mem_mb=1000 ,
           walltime= "01:00"
        shell:
           """

           minigraph -xasm -t 10 --call {input.graph} {input.assembly} > {output}

           """

def get_ordered_breeds(assemb_list,ref_genome):
    ordered_breed=[]
    ordered_breed.append(ref_genome)
    for anim in assemb_list:
        if anim != ref_genome:
            ordered_breed.append(anim)
    return ordered_breed


ordered_breed=get_ordered_breeds(assemb_list,ref_genome)

rule combine_per_chromo_minigraph_sv:
        input:expand("graph/minigraph/call_sv/sv_{{chromo}}_{breed}.bed",breed=ordered_breed)
        output:"graph/minigraph/call_sv/sv_{chromo}_combined.bed"
        threads:10
        resources:
           mem_mb= 500,
           walltime= "00:10"
        params:
            ref_bed="graph/minigraph/call_sv/sv_{chromo}_" + ref_genome + ".bed"
        shell:
           """

           #1       137187  138168  >s39    >s41    >s40:981:+:1:137180:138188

          paste <(awk '{{ print $1,$2,$3,$4,$5 }}' {params.ref_bed}) \
          <(awk '{{ L[FNR]=L[FNR] $NF "\\t" }} END {{ for(i=1;i<=FNR;i++)print L[i] }}' {input}) >> {output} 

           """

rule combine_all_sv_minigraph:
        input:expand("graph/minigraph/call_sv/sv_{chromo}_combined.bed",chromo=chromo_list)
        output:"graph/minigraph/call_sv/sv_all_call.bed"
        threads:2
        resources:
           mem_mb=2000 ,
           walltime= "01:00"
        params:
            breed_list=ordered_breed
        shell:
           """

           # add header 
           echo "chromo start_pos stop_pos start_node stop_node {params.breed_list}" >> {output}
           cat {input} >> {output}
           """


rule split_breakpoint_minigraph:
    input: "graph/minigraph/call_sv/sv_all_call.bed"
    output: 
        bed_left="graph/minigraph/call_sv/sv_left.bed",
        bed_right="graph/minigraph/call_sv/sv_right.bed"
    threads:10
    resources:
        mem_mb=2000,
        walltime= "00:10"
    shell:
        """
        sed 's/_UCD//g' {input}|
        awk 'NR>1{{cs++;print $1,$2,$2+2,"b"cs}}' OFS="\t" > {output.bed_left}


        sed 's/_UCD//g' {input}|
        awk 'NR>1{{cs++;print $1,$3,$3+2,"b"cs}}' OFS="\t" > {output.bed_right}


        """

rule annotate_sv_minigraph:
    input:"graph/minigraph/call_sv/sv_{side}.bed"
    output:
        raw="graph/minigraph/call_sv/sv_annot_raw_{side}.bed",
        summary="graph/minigraph/call_sv/sv_annot_summary_{side}.bed",
    threads: 10
    resources:
        mem_mb=2000,
        walltime="00:30"
    params:
        gff=GFF_file
    shell:
        """
        
        bedtools intersect -wb -loj -a {input} -b {params.gff} > {output.raw}

        ./annot_bubble.py -i {output.raw} -o {output.summary}

        """

rule combine_annotation_minigraph:
    input:expand("graph/minigraph/call_sv/sv_annot_summary_{side}.bed", side=["left","right"])
    output:"graph/minigraph/call_sv/sv_annot_combined.bed"
    threads: 1
    resources:
        mem_mb=1000,
        walltime="00:30"
    shell:
        """
        
        paste -- {input}|
        awk '{{print $1,$2,$3,$4,$5,$6,$9,$10,$11,$12}}'|
        sed '1i\svid chromo pos_left bp_left id_left name_left bp_right id_right name_right' |
        awk 'NR != 2'  > {output}

        """
