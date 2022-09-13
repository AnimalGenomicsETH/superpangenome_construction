
rule all:
    input:
        expand('VNTR_LR_genotyping/{sample}.bed',sample=config['samples'])

rule advntr_model:
    input:
        'VNTR_LR_genotyping/testable_VNTRs.txt'
    output:
        #'VNTR_LR_genotyping/test/models.db'
        'VNTR_LR_genotyping/d0_m100_L100/models.db'
    conda:
        'VNTR'
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '120:00'
    shell:
        '''
        while read -r -a p
        do
          advntr addmodel -r {config[reference]} -m {output} \
          -c "${p[0]}" -s "${p[1]}" -e "${p[2]}" -p "${p[3]}" -g "${p[4]}" \
          2> /dev/null
        done < {input}
        '''
rule advntr_genotype:
    input:
        bam = lambda wildcards: config['samples'][wildcards.sample],
        models = rules.advntr_model.output
    output:
        'VNTR_LR_genotyping/{sample}.bed'
    conda:
        'VNTR'
    threads: 2
    resources:
        mem_mb = 30000,
        walltime = '24:00',
        disk_scratch = 10
    shell:
        '''
        advntr genotype \
        -a {input.bam} \
        -m {input.models} \
        -p --naive\
        --working_directory $TMPDIR \
        -t {threads} \
        -of bed \
        -o {output} \
        --haploid
        '''
 
#manual paste to form csv 
#paste <(awk '{print $1"_"$2"\t"$7}' BSW.bed ) <(cut -f 8 BSW.bed) <(cut -f 8 Gaur.bed) <(cut -f 8 OxO1.bed) <(cut -f 8 OxO2.bed) <(cut -f 8 Pied.bed) <(cut -f 8 Nellore.bed) > truth_genotypes.csv