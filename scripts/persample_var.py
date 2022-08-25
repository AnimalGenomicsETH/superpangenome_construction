#!/usr/bin/env bash 

prog_list=["cactus","pggb"]
chromo_all=range(1,30)
chromo_list=[x for x in chromo_all if x not in [5,6,12,17,20,21]]
breed_list="OBV,BSW,Pied,Highland,Angus,Simmental,Brahman,Nellore,Gaur,Bison,Yak"
breed_list=breed_list.split(",")
assembly_dir="/cluster/work/pausch/to_share/sv_analysis/snarl/merged_assembly/vcf_split"

rule all:
    input:
        # expand("snarl/snarl_new/persample/{prog}/{prog}_{breed}_{chromo}.vcf",prog=prog_list,chromo=chromo_list,breed=breed_list),
        "snarl/snarl_new/persample/sv_persample_stat.tsv"


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

        bcftools view --min-ac 1 -I -a -s {wildcards.chromo}_{wildcards.breed} {input} |

        bcftools view -i "abs(ILEN)>=50" > {output}

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

rule count_breed_variants:
    input:
        expand("snarl/snarl_new/persample/{prog}/{prog}_{breed}_{chromo}.vcf",chromo=chromo_list,breed=breed_list,prog=prog_list),
        expand("snarl/snarl_new/persample/assembly/assembly_{breed}_{chromo}.vcf",breed=breed_list,chromo=chromo_list)
    output:
        "snarl/snarl_new/persample/sv_persample_stat.tsv"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        for file in {input}
        do
            prog=$(basename $file | cut -f1 -d"_")
            chromo=$(basename $file | cut -f2 -d"_")
            breed=$(basename $file | cut -f3 -d"_")
            echo $prog $chromo $breed $(wc -l < $file) >> {output}
        done

        """





