#!/usr/bin/env python


prog_list=["cactus","pggb"]
chromo_list=range(1,30)
snarldir="/cluster/work/pausch/to_share/sv_analysis/final_data/snarl/snarl_new"

rule all:
    input:
        # expand("normvar/graph/{prog}_{chromo}_norm.vcf",prog=prog_list,chromo=chromo_list),
        # expand("normvar/sr/var{chromo}_norm.vcf",chromo=chromo_list)
        # expand( "normvar/overlap/{prog}/{prog}_{chromo}_overlap_excl.tsv",chromo=chromo_list,prog=prog_list),
        # expand("normvar/overlap/{prog}/{prog}_{chromo}_overlap_nosvrep.tsv",chromo=chromo_list,prog=prog_list),
        "normvar/overlap/stat_fil.tsv"

rule norm_fill:
    input:
        snarldir + "/{prog}/{prog}_{chromo}_snarl_edit.vcf"
         # "{prog}{chromo}_var_graph.vcf"
    output:
        header="normvar/graph/{prog}_{chromo}_header.tsv",
        corvcf="normvar/graph/{prog}_{chromo}_cor.vcf",
        norm="normvar/graph/{prog}_{chromo}_norm.vcf"
    threads: 10
    resources:
        mem_mb= 2000 ,
        walltime= "04:00"
    shell:
        """

        awk '$1 ~/#/' {input} |
        awk '$0 !~ /ID=AT/{{print $0;next}}{{printf "%s\\n##INFO=<ID=CONFLICT,Number=1,Type=String,Description=Uncertain Paths>\\n",
$0}}' | sed 's/_UCD//g' > {output.header}

        bcftools reheader -h {output.header} <( sed 's/_UCD//g' {input}) > {output.corvcf}

        #bcftools view -i "abs(ILEN)<50" {output.corvcf} |
        #bcftools norm -f $REFGEN -m -any > {output.norm}

        awk '$1 ~ /#/ || (length($4) < 50 && length($5) < 50)' {output.corvcf} |
        bcftools norm -f $REFGEN -m -any > {output.norm}
        """


rule normalize_deepvar:
    input:
        "long_read_call/var{chromo}_pass.vcf.gz"
    output:
        "normvar/sr/var{chromo}_norm.vcf"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "04:00"
    shell:
        """

        #bcftools view -i "abs(ILEN)<50" {input} |

        awk '$1 ~ /#/ || (length($4) < 50 && length($5) < 50)'  {input}|
        bcftools norm -f $REFGEN -m -any > {output} 

        """

rule overlap_graph_sr:
    input:
        graph="normvar/graph/{prog}_{chromo}_norm.vcf",
        lr="normvar/sr/var{chromo}_norm.vcf"
    output:
        "normvar/overlap/{prog}/{prog}_{chromo}_overlap.tsv"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "02:00"
    shell:
        """
        
        ./concor_decons_vcf2.py {input.graph} BSW {input.lr} > {output}

        """


rule exclude_deepvar_near:
    input:
        "normvar/overlap/{prog}/{prog}_{chromo}_overlap.tsv"
    output:
        "normvar/overlap/{prog}/{prog}_{chromo}_overlap_excl.tsv"
    threads: 10
    resources:
        mem_mb= 500,
        walltime= "01:00"
    run:
        prev_pos=0
        with open(input[0]) as infile, open(output[0],"w") as outfile:
            for line in infile:
                token=line.strip().split()
                if abs(int(token[0])-prev_pos) < 50 and token[-1] == "no_match":
                    print(*token,"filtered",file=outfile)
                else:
                    print(*token,"keep",file=outfile)
                prev_pos=int(token[0])


rule exclude_assembly_sv:
    input:
        vcf="normvar/overlap/{prog}/{prog}_{chromo}_overlap_excl.tsv",
        svloc="svcoord_all_fixed.tsv"
    output:
        bedvar="normvar/overlap/{prog}/{prog}_{chromo}_overlap_nosv.bed",
        filter_var="normvar/overlap/{prog}/{prog}_{chromo}_overlap_nosv.tsv"
    threads: 1
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        awk -v chromo={wildcards.chromo} '{{print chromo,$1,$1+length($2),$2,$3,$(NF-1),$NF}}' OFS="\\t" {input.vcf} > {output.bedvar}


        bedtools subtract -A -a {output.bedvar} -b {input.svloc} > {output.filter_var}



        """


rule exclude_deep_variant_norep:
    input:
        vcf="normvar/overlap/{prog}/{prog}_{chromo}_overlap_nosv.tsv",
        svloc="ARS.repeats.sorted.bed"
    output:
        filter_var="normvar/overlap/{prog}/{prog}_{chromo}_overlap_nosvrep.tsv"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        bedtools subtract -A -a {input.vcf} -b {input.svloc} > {output.filter_var}



        """

localrules: combine_stat
rule combine_stat:
    input:
        expand("normvar/overlap/{prog}/{prog}_{chromo}_overlap_excl.tsv",prog=prog_list,chromo=chromo_list)
    output:
        "normvar/overlap/stat_fil.tsv"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "04:00"
    shell:
        """

        for file in {input}
        do 
            prog=$(basename $file | cut -f1 -d"_")
            chromo=$(basename $file | cut -f2 -d"_")

            total_var=$(wc -l < $file)
            total_match=$(grep -w "match" $file| wc -l)
            
            # SNP

            total_snp=$(awk 'length($2)==1 && length($3)==1' $file| wc -l )
            total_snp_match=$(grep -w "match" $file| awk 'length($2)==1 && length($3)==1'| wc -l)

            # Indel 
            total_indel=$(awk 'length($2)>1 || length($3)>1'  $file | wc -l)
            total_indel_match=$(awk 'length($2)>1 || length($3)>1'  $file | grep -w "match" | wc -l)
            
            # filtered 
            fil=$(grep -v "filtered" $file | wc -l)
            fil_match=$(grep -v "filtered" $file | grep -w "match"|wc -l)

            #SNP 

            fil_snp=$(awk 'length($2)==1 && length($3)==1' $file| grep -v "filtered" | wc -l )
            fil_snp_match=$(grep -w "match" $file| awk 'length($2)==1 && length($3)==1' | grep -v "filtered" | wc -l)
            
            #INDEL

            fil_indel=$(awk 'length($2)>1 || length($3)>1'  $file | grep -v "filtered" | wc -l)
            fil_indel_match=$(awk 'length($2)>1 || length($3)>1'  $file | grep -w "match" | grep -v "filtered" | wc -l)

            echo $prog $chromo $total_var $total_match $fil $fil_match \
                 $total_snp $total_snp_match $total_indel $total_indel_match \
                 $fil_snp $fil_snp_match $fil_indel $fil_indel_match >> {output}


        done 




        """

