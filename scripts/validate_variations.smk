breeds_list= ['Angus', 'Highland', 'OBV', 'Brahman', 'Yak', 'BSW', 'Pied', 'Gaur', 'Nellore', 'Simmental', 'Bison']
fasta_dir="/cluster/work/pausch/danang/psd/scratch/real_comp/assembly"
snarl_dir="/cluster/work/pausch/danang/psd/scratch/real_comp/graph/snarl"
chromo_list=range(1,30)
prog_list=["cactus","pggb"]

wildcard_constraints:
    breed = r'|'.join(breeds_list)

rule all:
    input:
        expand('assembly_{chromo}.wave.SV.vcf',chromo=chromo_list)
        # expand("{breeds}/{breeds}_overlap_stat",breeds=breeds_list),
        #expand("{breeds}/{breeds}_{chromo}_{prog}_paf_overlap_detail.tsv",breeds=breeds_list,chromo=chromo_list,prog=prog_list),
        #expand("{breeds}/{breeds}_overlap_snp_stat.tsv",breeds=breeds_list),
        #expand("comb_stat_snp_{prog}_overlap.tsv",prog=prog_list),
        #expand("comb_stat_snp_{prog}_overlap_filtered.tsv",prog=prog_list)

asm_map = {'Angus':'asm5','Highland':'asm5','OBV':'asm5','BSW':'asm5','Pied':'asm5','Simmental':'asm5','Nellore':'asm10','Brahman':'asm10','Gaur':'asm20','Bison':'asm20','Yak':'asm20'}

rule map_assembly:
    input:
        ref_genome= fasta_dir + "/{chromo}/UCD_{chromo}.fa",
        alt_genome = fasta_dir + "/{chromo}/{breed}_{chromo}.fa"
    output: "{breed}_{chromo}.paf"
    params:
        preset = lambda wildcards: asm_map[wildcards.breed]
    threads: 1
    resources:
        mem_mb= 10000,
        walltime= "02:00"
    shell:
        """
        minimap2 -cx {params.preset} -t{threads} \
        --cs {input.ref_genome} {input.alt_genome} |\
        sort -k6,6 -k8,8n > {output}
        """


rule call_variations:
    input:
        paf = "{breed}_{chromo}.paf",
        ref = fasta_dir + "/{chromo}/UCD_{chromo}.fa"
    output:
        "{breed}_{chromo}.vcf"
    threads: 1
    resources:
        mem_mb= 2000,
        walltime= "30"
    shell:
        """
        paftools.js call -f {input.ref} -s {wildcards.breed} {input.paf} > {output}
        """

#concat and then merge
rule bcftools_concat:
    input:
        expand("{{breed}}_{chromo}.vcf",chromo=range(1,30))
    output:
        "{breed}.vcf"
    shell:
        'bcftools concat -o {output} {input}'

rule bcftools_merge:
    input:
        expand("{breed}_{{chromo}}.vcf",breed=breeds_list)
    output:
        'assembly_{chromo}.vcf'
    params:
        ' '.join(expand("{breed}_{{chromo}}.vcf.gz",breed=breeds_list))
    shell:
        '''
        for i in {input}
        do
        bgzip -c $i > $i.gz
        tabix -fp vcf $i.gz
        done
        bcftools merge -o {output} {params}
        '''

rule assembly_SVs:
    input:
        vcf = 'assembly_{chromo}.vcf',
        ref = fasta_dir + "/{chromo}/UCD_{chromo}.fa"
    output:
        'assembly_{chromo}.wave.SV.vcf'
    resources:
        mem_mb = 2500
    shell:
        '''
        vcfwave {input.vcf} |\
        bcftools norm -m - -f {input.ref} |\
        bcftools view -o {output} -i 'abs(ILEN)>=50'
        '''
        
rule overlap_variations:
    input:
        snarl_var = snarl_dir + "/{prog}/{prog}_{chromo}_snarl.vcf",
        paf_var = "{breeds}/{breeds}_{chromo}_varassemb.tsv"
    output: "{breeds}/{breeds}_{chromo}_{prog}_paf_overlap.tsv"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "00:30"
    shell:
        """
        
         ./con_dec.py {input.snarl_var} {wildcards.breeds} {input.paf_var} |
         awk -v prog={wildcards.prog} -v chromo={wildcards.chromo} '{{print prog,chromo,$0}}' OFS="\\t" > {output}

        """

rule combine_stat:
    input:expand("{{breeds}}/{{breeds}}_{chromo}_{prog}_paf_overlap.tsv",chromo=chromo_list,prog=prog_list)
    output:"{breeds}/{breeds}_overlap_stat"
    threads: 10
    resources:
        mem_mb= 1000 ,
        walltime= "00:10"
    shell:
        """
        
        cat {input} > {output}

        """


rule overlap_details:
    input:
        snarl_var = snarl_dir + "/{prog}/{prog}_{chromo}_snarl.vcf",
        paf_var = "{breeds}/{breeds}_{chromo}_varassemb.tsv"
    output: "{breeds}/{breeds}_{chromo}_{prog}_paf_overlap_detail.tsv"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "00:30"
    shell:
        """
        
         ./con_details.py {input.snarl_var} {wildcards.breeds} {input.paf_var} > {output}

        """

rule overlap_snps:
    input:
        snarl_var = snarl_dir + "/{prog}/{prog}_{chromo}_snarl.vcf",
        paf_var = "{breeds}/{breeds}_{chromo}_varassemb.tsv"
    output: "{breeds}/{breeds}_{chromo}_{prog}_paf_overlap_snp.tsv"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "00:30"
    shell:
        """
        
         ./con_dec.py {input.snarl_var} {wildcards.breeds} {input.paf_var} |
         awk -v prog={wildcards.prog} -v chromo={wildcards.chromo} '{{print prog,chromo,$0}}' OFS="\\t" > {output}

        """

rule combine_stat_snps:
    input:expand("{{breeds}}/{{breeds}}_{chromo}_{prog}_paf_overlap.tsv",chromo=chromo_list,prog=prog_list)
    output:"{breeds}/{breeds}_overlap_snp_stat.tsv"
    threads: 10
    resources:
        mem_mb= 1000 ,
        walltime= "00:10"
    shell:
        """
        
        cat {input} > {output}

        """


rule combine_stat_details:
    input:
        expand("{breeds}/{breeds}_{chromo}_{{prog}}_paf_overlap_detail.tsv",breeds=breeds_list,chromo=chromo_list,prog=prog_list)
    output:
        "comb_stat_snp_{prog}_overlap.tsv"
    threads: 10
    resources:
        mem_mb= 1000 ,
        walltime= "04:00"
    shell:
        """

        for file in {input}
        do
            breed=$(cut -f1 -d"_" <(basename $file))
            chromo=$(cut -f2 -d"_" <(basename $file))
            prog=$(cut -f3 -d"_" <(basename $file))
            awk -v breed=$breed -v chromo=$chromo -v prog=$prog \
                    '{{ print breed,chromo,prog,$0 }}' OFS="\\t" $file >> {output}
        done

        """


rule exclude_near:
    input:
        "comb_stat_snp_{prog}_overlap.tsv"
    output:
        "comb_stat_snp_{prog}_overlap_filtered.tsv"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    run:
        prev_pos=0
        with open(input[0]) as infile, open(output[0],"w") as outfile:
            for line in infile:
                token=line.strip().split()
                if abs(int(token[3])-prev_pos) < 50 and token[6] == "no_match":
                    print(*token,"filtered",file=outfile)
                else:
                    print(*token,"keep",file=outfile)
                prev_pos=int(token[3])

