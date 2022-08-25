#!/usr/bin/env python


chromo_all=list(range(1,30))
# chromo_list.append(19)
prog_list=["cactus","pggb"]
chromo_rem=[5,6,12,17,20,21]
#chromo_list=[12,2,4,20,16,17,18,21,5,6]
chromo_list=[x for x in chromo_all if x not in chromo_rem]
breed_list="OBV,BSW,Pied,Highland,Angus,Simmental,Brahman,Nellore,Gaur,Bison,Yak"
breed_list=breed_list.split(",")
assembly_dir="/cluster/work/pausch/to_share/sv_analysis/snarl/merged_assembly/vcf_split"

intersect_param={
        "A":"max_dist=25",
        "B":"max_dist_linear=0.5 max_dist=1000",
        "C":"max_dist_linear=0.5 min_seq_id=0.5"
        }


rule all:
    input:
  #      expand("paramtr/{chromo}_{grpar}_merged.vcf",chromo=chromo_list,grpar=intersect_param.keys()),
          # expand("paramtr/{grpar}_stat_merged.vcf",grpar=intersect_param.keys()),
          # expand("snarl/snarl_new/persample/{prog}/{prog}_{breed}_{chromo}.vcf",prog=prog_list,chromo=chromo_list,breed=breed_list),
          # expand("snarl/snarl_new/persample/assembly/assembly_{breed}_{chromo}.vcf",breed=breed_list,chromo=chromo_list)
          # expand("paramtr/persample/{breed}/{breed}_{chromo}_merged.vcf",breed=breed_list,chromo=chromo_list)
          "paramtr/breed_stat_merged_breed.tsv"


localrules: create_overlap_file
rule create_overlap_file:
    input:
        expand("snarl/snarl_new/{prog}/{prog}_{{chromo}}_snarl_norm.vcf",prog=prog_list),
        "/cluster/work/pausch/to_share/sv_analysis/snarl/merged_assembly/vcf_split/{chromo}_assembly_norm_fil.vcf"
    output:
        listfiles="paramtr/{chromo}_merged.tsv"
    shell:
        """

        echo {input} | tr ' ' '\\n' > {output.listfiles}

        """

rule merge_jasmine:
    input:
        "paramtr/{chromo}_merged.tsv"
    output:
        "paramtr/{chromo}_{grpar}_merged.vcf"
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
        expand("paramtr/{chromo}_{{grpar}}_merged.vcf",chromo=chromo_list)
    output:
        "paramtr/{grpar}_stat_merged.vcf"
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


rule persample_var:
    input:
        "snarl/snarl_new/{prog}/{prog}_{chromo}_snarl_norm.vcf"
    output:
        "snarl/snarl_new/persample/{prog}/{prog}_{breed}_{chromo}.vcf"
    threads: 5
    resources:
        mem_mb= 1000,
        walltime= "01:00"
    shell:
        """

        bcftools view --min-ac 1 -I -a -s {wildcards.chromo}_{wildcards.breed} {input}|
        bcftools view -i "abs(ILEN)>=50"> {output}

        """

rule assembly_breed_var:
    input:
        assembly_dir + "/{chromo}_assembly_norm_fil.vcf"
    output:
        "snarl/snarl_new/persample/assembly/assembly_{breed}_{chromo}.vcf"
    threads: 2
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    shell:
        """

        bcftools view --min-ac 1 -I -a -s {wildcards.breed} {input} | 
        bcftools view -i "abs(ILEN)>=50" > {output}
        

        """


localrules: create_overlap_brfile
rule create_overlap_brfile:
    input:
        expand("snarl/snarl_new/persample/{prog}/{prog}_{{breed}}_{{chromo}}.vcf",prog=prog_list),
        "snarl/snarl_new/persample/assembly/assembly_{breed}_{chromo}.vcf"
    output:
        listfiles="paramtr/persample/{breed}/{breed}_{chromo}_mergedlist.tsv"
    shell:
        """

        echo {input} | tr ' ' '\\n' > {output.listfiles}

        """

rule merge_jasmine_breed:
    input:
        "paramtr/persample/{breed}/{breed}_{chromo}_mergedlist.tsv"
    output:
        "paramtr/persample/{breed}/{breed}_{chromo}_merged.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        tempfil=/cluster/scratch/cdanang/jastemp/jastemp_{wildcards.chromo}_{wildcards.breed}

        if [[ ! -d $tempfil ]]; then
            mkdir $tempfil
        fi

        jasmine file_list={input} out_dir=$tempfil \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_reads=0 max_dist_linear=0.5 max_dist=1000  \
        --ignore_strand --threads=8 --allow_intrasample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output}


        """


rule combine_information_breed:
    input:
        expand("paramtr/persample/{breed}/{breed}_{chromo}_merged.vcf",chromo=chromo_list,breed=breed_list)
    output:
        "paramtr/breed_stat_merged_breed.tsv"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

         for file in {input}
         do
            breed=$(basename $file | cut -f1 -d"_")
            chromo=$(basename $file | cut -f2 -d"_")
            grep -v "#" $file |
            awk -v breed=$breed -v chromo=$chromo '{{match($0,/Name=([^;]+)/,arr);match($0,/SUPP_VEC=([0-9]+)/,vec);print $1,$2,vec[1],arr[1],breed,chromo}}' OFS="\t" >> {output}
          done


        """
