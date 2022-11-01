rule TRF:
    input:
        'assemblies/{chromosome}/{sample}.fa'
    output:
        temp('assemblies/{chromosome}/{sample}.fa.2.7.7.80.10.50.500.dat')
    threads: 1
    resources:
        mem_mb = 1000,
        walltime = '4:00'
    params:
        lambda wildcards, output: ' '.join(output[0].split('.')[-8:-1])
    shell:
        '''
        TRF {input} {params} -h
        '''

rule postprocess_TRF:
    input:
        TR = rules.TRF.output[0],
        repeats = 'assemblies/{chromosome}/{sample}.fa_rm.bed'
    output:
        'VNTRs/TRF/{sample}.{chromosome}.TR.bed'
    params:
        min_length = 100,
        max_length = 10000,
        min_TR_length = 10,
        max_TR_length = 100
    shell:
        '''
        awk -v m={params.min_length} -v M={params.max_length} -v l={params.min_TR_length} -v L={params.max_TR_length} '{{if($1~/Sequence/){{C=$2}} else{{ if($1~/[[:digit:]+]/&&($2-$1)>m&&($2-$1)<M&&length($14)>=l&&length($14)<=L) {{print C"\\t"$1"\\t"$2"\\t"$8"\\t"$14}}}}}}' {input} |\
        grep -P "^\d" |\
        sort -k1,1n -k2,2n |\
        bedtools merge -c 4,5 -o collapse -i - |\
        bedtools subtract -A -f 0.9 -F 0.000001 -a {input.TRs} -b {input.repeats} > {output.TR} |\
        python postprocess_TRs.py > {output}
        '''


rule odgi_position:
    input:
        gfa = 'graphs/{pangenome}/{chromosome}.gfa',
        TRs = expand('VNTRs/TRF/{sample}.{{chromosome}}.TRs.bed',sample=get_reference_ID())
    output:
        'VNTRs/{pangenome}/{chromosome}.liftover.bed'
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '4:00'
    shell:
        '''        
        odgi position -i {input.gfa} -b {input.TRs} -t {threads} > {output}
        '''

## simplify bubble -> coords
rule gfatools_bubble:
    input:
        'graphs/minigraph/{chromosome}.gfa'
    output:
        'VNTRs/minigraph/{chromosome}.bubbles.bed'
    shell:
        '''
        gfatools bubble {input} | cut -f -3,12 > {output}
        '''

rule estimate_minigraph_coordinates:
    input:
        bubbles = rules.gfatools_bubble.output[0],
        TRs = expand('VNTRs/TRF/{sample}.{{chromosome}}.TRs.bed',sample=get_reference_ID())
    output:
        'VNTRs/minigraph/{chromosome}.liftover.bed'
    params:
        threshold = 0.000000001, #lambda wildcards: int(wildcards.rate)/100,
        dist = config.get('merging_distance',0)
    shell:
        '''
        bedtools intersect -f {params.threshold} -wo -a {input.bubbles} -b {input.TRs} |\
        awk '{{n=split($4, a, ","); print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t"$4"\\t"n"\\t"$8}}' |\
        bedtools merge -d {params.dist} -c 4,5,6,7,8 -o min,max,distinct,sum,distinct -delim '|'  -i {input} > {output}
        '''

import regex
from collections import defaultdict
def count_VNTRs(sequences,TR):
    allowed_errors = int(len(TR)*config.get('TR_divergence_limit',0)) if config.get('scaling','linear') != 'log' else math.floor(math.sqrt(len(TR)))
    counts = dict()
    for asm, sequence in sequences.items():
        asm_ID = asm.split('_')[1].split(':')[0] if '_' in asm else asm
        if len(sequence) > config.get('max_region_length',100000):
            counts[asm_ID] = -2
        try:
            counts[asm_ID] = len(regex.findall(f"(?eV1)({TR}){{2i+2d+1s<={allowed_errors}}}",sequence,concurrent=True,timeout=config.get('timeout',60)))
        except TimeoutError:
            counts[asm_ID] = -1
            break
    return counts


rule process_VNTRs:
    input:
        ''
    output:
        'bad'
    run:
        pass


rule advntr_model:
    input:
        'VNTR_LR_genotyping/testable_VNTRs.txt'
    output:
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
        scratch = '10G'
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
