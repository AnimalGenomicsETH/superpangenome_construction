!/usr/bin/env python

prog_list=["cactus","pggb","minigraph"]
#chromo_list=list(range(1,4)) +  list(range(5,30)) 
chromo_list=list(range(1,30))
for ch in [4,5,6,21]:
    chromo_list.remove(ch)

rule all:
    input:
        expand("{prog}/{prog}_{chromo}_sv_filter.vcf",prog=prog_list,chromo=chromo_list),
        expand("{prog}/{prog}_{chromo}_sv_interval.vcf",prog=prog_list,chromo=chromo_list),
        expand("merged/{chromo}_merged_all.vcf",chromo=chromo_list),
        "merged/merged_overlap.tsv",
        expand("merged/{chromo}_sv_interval_overlap.vcf",chromo=chromo_list),
        expand("merged/{chromo}_merged_intjas.vcf",chromo=chromo_list),
        "merged/merged_overlap_intjas.tsv",
        expand("{prog}/{prog}_{chromo}_sv_wave.vcf",prog=prog_list,chromo=chromo_list),
        expand( "{prog}/{prog}_{chromo}_sv_wavefil.vcf",prog=prog_list,chromo=chromo_list),
        expand("merged/{chromo}_merged_wave.vcf",chromo=chromo_list),
        "merged/merged_overlap_wave.tsv",
        expand("{prog}/{prog}_{chromo}_sv_wavenorep.vcf",chromo=chromo_list,prog=prog_list),
        expand("merged/{chromo}_merged_wave_norep.vcf",chromo=chromo_list),
        expand("merged_assembly/{chromo}_assembly_wave_merged.vcf",chromo=chromo_list),
        "merged_assembly/merged_overlap_assembly_wave.tsv",
        expand("merged/annot/{chromo}_merged_wave_annot.vcf",chromo=chromo_list),
        "merged/annot/annot_merged_norep_tidy.tsv"



rule filter_vcf:
    input:
        "{prog}/{prog}_{chromo}_sv_trim2.vcf"
    output:
        "{prog}/{prog}_{chromo}_sv_filter.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        awk '$1 ~ /^#/ || length($4) <= 100000 || length($5) <= 100000' {input} > {output}
        """


rule merge_jasmine:
    input:
        expand("{prog}/{prog}_{{chromo}}_sv_filter.vcf",prog=prog_list)
    output:
        listfiles="merged/{chromo}_all_list.tsv",
        merged="merged/{chromo}_merged_all.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        echo {input} | tr ' ' '\n' > {output.listfiles}


        jasmine file_list={output.listfiles} \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_reads=0 \
        --ignore_strand --threads=8 --allow_intrasample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output.merged}

        """


rule combine_overlap:
    input:
        expand("merged/{chromo}_merged_all.vcf",chromo=chromo_list)
    output:
         "merged/merged_overlap.tsv"
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

rule create_interval:
    input:
        "{prog}/{prog}_{chromo}_sv_filter.vcf"
    output:
        "{prog}/{prog}_{chromo}_sv_interval.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """
         grep -v "#" {input}  |
         awk '{{ print $1,$2,$2+length($4),$1"_"$2"_"$2+length($4) }}' OFS="\\t" > {output}
        """


rule merge_interval:
    input:
        expand("{prog}/{prog}_{{chromo}}_sv_interval.vcf",prog=prog_list)
    output:
        "merged/{chromo}_sv_interval_overlap.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        bedtools multiinter -header -i {input} -names cactus pggb minigraph > {output}
        #bedtools intersect -a {input[0]} -b {input[1]} {input[2]} \
        #        -sorted -wa -wb -names pggb minigraph > {output}

        """


rule merge_jasmine_interval:
    input:
        listfiles="merged/{chromo}_all_list.tsv",
        inputfiles=expand("{prog}/{prog}_{{chromo}}_sv_filter.vcf",prog=prog_list)
    output:
        merged="merged/{chromo}_merged_intjas.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        jasmine file_list={input.listfiles} \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_reads=0 max_dist=25 --use_end \
        --ignore_strand --threads=8 --allow_intrasample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output.merged}

        """


rule combine_overlap_intjas:
    input:
        expand("merged/{chromo}_merged_intjas.vcf",chromo=chromo_list)
    output:
        "merged/merged_overlap_intjas.tsv"
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


rule vcf_wave:
    input:
         "{prog}/{prog}_{chromo}_sv_filter.vcf"
    output:
         "{prog}/{prog}_{chromo}_sv_wave.vcf"
    threads: 10
    resources:
        mem_mb= 10000,
        walltime= "24:00"
    shell:
        """
        vcfwave -I 1000 -j 10 {input} > {output}
        """

rule filter_sv:
    input:
        "{prog}/{prog}_{chromo}_sv_wave.vcf"
    output:
        "{prog}/{prog}_{chromo}_sv_wavefil.vcf"
    threads: 2
    resources:
        mem_mb= 2000 ,
        walltime= "01:00"
    shell:
        """
         bcftools view -i 'abs(ILEN)>=50' {input} > {output}
        """

rule merge_wave:
    input:
        expand("{prog}/{prog}_{{chromo}}_sv_wavefil.vcf",prog=prog_list)
    output:
        listfiles="merged/{chromo}_wave_list.tsv",
        merged="merged/{chromo}_merged_wave.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        echo {input} | tr ' ' '\n' > {output.listfiles}


        jasmine file_list={output.listfiles} \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_reads=0 max_dist=25 --use_end \
        --ignore_strand --threads=8 --allow_intrasample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output.merged}

        """


rule exclude_rep:
    input:
        "merged/{chromo}_merged_wave.vcf"
    output:
        "merged/{chromo}_merged_wave_norep.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """
        
        bcftools view -T ^ARS.repeats.sorted.bed {input} > {output}

        """


rule combine_overlap_wave:
    input:
        expand("merged/{chromo}_merged_wave.vcf",chromo=chromo_list)
    output:
        "merged/merged_overlap_wave.tsv"
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


## exclude rep region

rule exclude_wave_rep:
    input:
        "{prog}/{prog}_{chromo}_sv_wavefil.vcf"
    output:
        "{prog}/{prog}_{chromo}_sv_wavenorep.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        bcftools view -T ^ARS.repeats.sorted.bed {input} > {output}

        """


rule merge_assembly_wave:
    input:
        expand("{prog}/{prog}_{{chromo}}_sv_wavefil.vcf",prog=prog_list),
        "merged_assembly/vcf_split/{chromo}_assembly.vcf"
    output:
        listfiles="merged_assembly/vcf_split/{chromo}_assembly_wave_merged.tsv",
        merged="merged_assembly/{chromo}_assembly_wave_merged.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        echo {input} | tr ' ' '\\n' > {output.listfiles}


        jasmine file_list={output.listfiles} \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_reads=0 max_dist=25 --use_end \
        --ignore_strand --threads=8 --allow_intrasample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output.merged}


        """


rule combine_overlap_assembly_wave:
    input:
        expand("merged_assembly/{chromo}_assembly_wave_merged.vcf",chromo=chromo_list)
    output:
        "merged_assembly/merged_overlap_assembly_wave.tsv"
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

rule annot_merge:
    input:
        "merged/{chromo}_merged_wave_norep.vcf"
    output:
        "merged/annot/{chromo}_merged_wave_nonrep_annot.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        singularity exec \
        --bind $PWD:/data \
        --bind /cluster/work/pausch/inputs/ref/BTA/UCD1.2:/ref \
        --bind /cluster/work/pausch/danang/psd/scratch/real_comp/graph/snarl:/annot \
        $SIFDIR/vep.sif \
        /opt/vep/src/ensembl-vep/vep \
        --force_overwrite --vcf --overlaps \
        --no_check_variants_order  \
        -i /data/{input} -o /data/{output} \
        -gff /annot/bos_annot.gff.gz -fasta /ref/ARS-UCD1.2_Btau5.0.1Y.fa

        """



rule annot_tidy:
    input:
        expand("merged/annot/{chromo}_merged_wave_nonrep_annot.vcf",chromo=chromo_list)
    output:
        "merged/annot/annot_merged_norep_tidy.tsv"
    threads: 10
    resources:
        mem_mb= 500,
        walltime= "01:00"
    shell:
        """
        
        for input_file in {input}
        do 
             extract_annot_full.py $input_file >> {output}
        done
        """


rule normalized_assembly:
    input:
         "{prog}/{prog}_{chromo}_sv_wave.vcf"
    output:
         "{prog}/{prog}_{chromo}_sv_wave_norm.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "24:00"
    shell:
        """

        vcfwave -I 1000 -j 10 {input} > {output}

        """



