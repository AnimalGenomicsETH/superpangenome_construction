#!/usr/bin/env python


chromo_all=list(range(1,30))
prog_list=["cactus","pggb","minigraph"]
chromo_rem=[5,6,12,17,20,21]
#chromo_list=[12,2,4,20,16,17,18,21,5,6]
chromo_list=[x for x in chromo_all if x not in chromo_rem]


intersect_param={
        "A":"max_dist=25",
        "B":"max_dist_linear=0.5 max_dist=1000",
        "C":"max_dist_linear=0.5 min_seq_id=0.5"
        }


rule all:
    input:
        expand("snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_norm.vcf",prog=prog_list,chromo=chromo_list)
        #expand("paramsearch/{grpar}_stat_merged.vcf",grpar=intersect_param.keys())
        # expand("graph/cactus/cactus_{chromo}.gfa",chromo=chromo_list)
        # expand("snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_edit.vcf",chromo=chromo_list,prog=prog_list)
 #       expand("paramsearch/{chromo}_merged.tsv",chromo=chromo_list),
  #      expand("paramsearch/{chromo}_{grpar}_merged.vcf",chromo=chromo_list,grpar=intersect_param.keys()),
 	 #expand("paramsearch/{grpar}_stat_merged.vcf",grpar=intersect_param.keys())


rule correct_cactus_path:
    input:
        "graph/cactus/cactus_simple_{chromo}.gfa"
    output:
        "graph/cactus/cactus_{chromo}.gfa"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        awk '{{if($1~/P/){{split($2,a,"\.");$2=a[2];print $0}}else{{print $0}}}}' OFS="\\t" {input} > {output}


        """


rule deconstruct_vcf:
    input:
        "graph/{prog}/{prog}_{chromo}.gfa"
    output:
        "snarl/snarl_new/{prog}/{prog}_{chromo}_snarl.vcf"
    threads: 10
    resources:
        mem_mb= 5000,
        walltime= "04:00"
    shell:
        """

        if [[ {wildcards.prog} = "minigraph"  ]]
        then 

        vg deconstruct -p {wildcards.chromo}_UCD -a -t {threads} \
        -v {input} > {output}

        else 

        vg deconstruct -p {wildcards.chromo}_UCD -a -e -d 1 -t {threads} \
        -v {input} > {output}
        
        fi

        """


rule edit_sed:
    input:
        "snarl/snarl_new/{prog}/{prog}_{chromo}_snarl.vcf"
    output:
        "snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_edit.vcf"
    threads: 2
    resources:
        mem_mb= 1000,
        walltime= "01:00"
    shell:
        """
        
        # sed 's/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g' {input} | bcftools view -a -i "abs(ILEN) < 100000"  > {output}


        if [[ {wildcards.prog} = "minigraph"  ]]
        then 

        sed 's/FORMAT//g' {input} > {output}

        else 

        sed 's/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t$/\\t./g' {input}  |\
        awk '$0 !~ /ID=AT/{{print $0;next}}{{printf "%s\\n##INFO=<ID=CONFLICT,Number=1,Type=String,Description=Uncertain Paths>\\n", $0}}' |\
         awk 'NF==20 || $1 ~ /#/' | bcftools view -a -i "abs(ILEN)<100000"  > {output}

        fi
        """



rule vcf_wave:
    input:
        "snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_edit.vcf"
    output:
        "snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_wave.vcf"
    threads: 10
    resources:
        mem_mb= 5000 ,
        walltime= "120:00"
    shell:
        """

         vcfwave -I 1000 -j 10 {input} > {output}

        """



# rule norm_fill:
    # input:
        # "snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_wave.vcf"
    # output:
        # header="snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_header.tsv",
        # corvcf="snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_cor.vcf",
        # norm="snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_norm.vcf"
    # threads: 10
    # resources:
        # mem_mb= 2000 ,
        # walltime= "04:00"
    # shell:
        # """

        # awk '$1 ~/#/' {input} |
        # awk '$0 !~ /ID=AT/{{print $0;next}}{{printf "%s\\n##INFO=<ID=CONFLICT,Number=1,Type=String,Description=Uncertain Paths>\\n", $0}}' > {output.header}

        # bcftools reheader -h {output.header} {input} > {output.corvcf}

        # bcftools norm -f $REFGEN -m -any {output.corvcf} |
        # bcftools view -i "abs(ILEN)>=50" > {output.norm}

        # """

rule norm_fill:
    input:
        "snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_wave.vcf"
    output:
        norm="snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_norm.vcf"
    threads: 10
    resources:
        mem_mb= 2000 ,
        walltime= "04:00"
    shell:
        """

        sed -E 's/CONFLICT=(.*);//g' {input}  | 
        sed 's/_UCD//g' | 
        bcftools norm -f $REFGEN -m -any |
        bcftools view -i "abs(ILEN)>=50" > {output.norm}

        """

localrules: create_overlap_file
rule create_overlap_file:
    input:
        expand("snarl/snarl_new/{prog}/{prog}_{{chromo}}_snarl_norm.vcf",prog=prog_list),
        "/cluster/work/pausch/to_share/sv_analysis/snarl/merged_assembly/vcf_split/{chromo}_assembly_norm_fil.vcf"
    output:
        listfiles="paramsearch/{chromo}_merged.tsv"
    shell:
        """

        echo {input} | tr ' ' '\\n' > {output.listfiles}

        """

rule merge_jasmine:
    input:
        "paramsearch/{chromo}_merged.tsv"
    output:
        "paramsearch/{chromo}_{grpar}_merged.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    params:
        gr=lambda w: intersect_param[w.grpar]
    shell:
        """

        tempfil=/cluster/scratch/cdanang/jastemp/jastemp_{wildcards.chromo}_{wildcards.grpar}

        if [[ ! -d $tempfil ]]; then
            mkdir $tempfil
        fi

        jasmine file_list={input} out_dir=$tempfil \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_reads=0 {params.gr} \
        --ignore_strand --threads=8 --allow_intrasample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output}


        """

rule combine_information:
    input:
        expand("paramsearch/{chromo}_{{grpar}}_merged.vcf",chromo=chromo_list)
    output:
        "paramsearch/{grpar}_stat_merged.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

         for file in {input}
         do
            grep -v "#" $file |
            awk '{{match($0,/Name=([^;]+)/,arr);match($0,/SUPP_VEC=([0-9]+)/,vec);print $1,$2,vec[1],arr[1]}}' OFS="\t" >> {output}
          done


        """


### per sample variations 

