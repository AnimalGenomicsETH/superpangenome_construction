rule repeatmasker_to_bed:
    input:
        'assemblies/{chromosome}/{sample}.fa.out'
    output:
        'assemblies/{chromosome}/{sample}.fa_rm.bed'
    shell:
        '''
        python /cluster/work/pausch/alex/software/RepeatMasker/util/RM2Bed.py {input} {output}
        '''

rule TRF:
    input:
        'assemblies/{chromosome}/{sample}.fa'
    output:
        temp('VNTRs/TRF/{sample}.{chromosome}.TRF')
    threads: 1
    resources:
        mem_mb = 1000,
        walltime = '4:00'
    shell:
        '''
        TRF {input} 2 7 7 80 10 50 500 -h -ngs > {output}
        '''
q
rule postprocess_TRF:
    input:
        TR = rules.TRF.output[0],
        repeats = 'assemblies/{chromosome}/{sample}.fa_rm.bed'
    output:
        TR = 'VNTRs/TRF/{sample}.{chromosome}.TR.bed'
    params:
        min_length = 100,
        max_length = 10000,
        min_TR_length = 10,
        max_TR_length = 100,
        TR_merge_window = 100
    shell:
        '''
        awk -v m={params.min_length} -v M={params.max_length} -v l={params.min_TR_length} -v L={params.max_TR_length} '{{if($1~/@/){{C=substr($1,2)}} else{{ if($1~/[[:digit:]+]/&&($2-$1)>m&&($2-$1)<M&&length($14)>=l&&length($14)<=L) {{print C"\\t"$1"\\t"$2"\\t"$8"\\t"$14}}}}}}' {input.TR} |\
        bedtools subtract -A -f 0.9 -F 0.000001 -a - -b {input.repeats} |\
        sort -k1,1 -k2,2n |\
        bedtools merge -d {params.TR_merge_window} -c 5 -o collapse |\
        awk -v OFS='\t' '{{split($4,st,","); min=""; for (i in st) {{  if (min == "" || length (st[i]) < length (min)) {{ min = st[i] }} }} print $1,$2,$3,min}}' > {output.TR}
        '''

rule odgi_position:
    input:
        gfa = 'graphs/{pangenome}/{chromosome}.gfa',
        TRs = expand('VNTRs/TRF/{sample}.{{chromosome}}.TR.bed',sample=get_reference_ID())
    output:
        'VNTRs/{pangenome,pggb|cactus}/{chromosome}.liftover.bed'
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
        'graphs/minigraph/{chromosome}.basic.gfa'
    output:
        'VNTRs/minigraph/{chromosome}.bubbles.bed'
    shell:
        '''
        gfatools bubble {input} | cut -f -3,12 > {output}
        '''

rule estimate_minigraph_coordinates:
    input:
        bubbles = rules.gfatools_bubble.output[0],
        TRs = expand('VNTRs/TRF/{sample}.{{chromosome}}.TR.bed',sample=get_reference_ID())
    output:
        'VNTRs/minigraph/{chromosome}.liftover.bed'
    params:
        threshold = 0.000000001, #lambda wildcards: int(wildcards.rate)/100,
        dist = config.get('merging_distance',0)
    shell:
        '''
        bedtools intersect -f {params.threshold} -wo -a {input.bubbles} -b {input.TRs} |\
        awk '{{n=split($4, a, ","); print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t"$4"\\t"n"\\t"$8}}' |\
        bedtools merge -d {params.dist} -c 4,5,6,7,8 -o min,max,distinct,sum,distinct -delim '|' > {output}
        '''


while read a; do awk -v X=$a -v SI=$(echo $a | cut -d',' -f 1) -v SO=$(echo $a | rev| cut -d',' -f 1 | rev) '$4==(">"SI)&&$5==(">"SO)&&$6!="." {split($6,a,":"); print a[4]":"a[5]"-"a[6]}' graphs/minigraph/*bed | tr '\n' '\t'; echo -en "\n" ; done < <(cut -f 6 VNTRs/minigraph/28.liftover.bed | head)


chrom, start, end, TR_start, TR_end, nodes, TR = line.rstrip().split()
    count = nodes.count(',')+1
    if '|' in TR:
        TR = Levenshtein.median(TR.split('|'))
    if config.get('low',2) < count < config.get('high',math.inf) and config.get('TR_low',5) <= len(TR) <= config.get('TR_high',math.inf):

        regions = []
        for node in nodes.split('|'):
            source, *_, sink = node.split(',')
            regions.extend(subprocess.run(f"awk '$4==\">{source}\"&&$5==\">{sink}\"' {chrom}/*bed",shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1])


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
