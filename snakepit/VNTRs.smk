pangenome_samples = config['pangenome_samples']
def get_reference_ID():
    return 'HER'

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
    threads: 2
    resources:
        mem_mb = 2500,
        walltime = '4:00'
    shell:
        '''        
        odgi position -i {input.gfa} -b {input.TRs} -t {threads} > {output}
        '''

## simplify bubble -> coords
rule gfatools_bubble:
    input:
        gfa = 'graphs/minigraph/{chromosome}.basic.gfa',
        TRs = expand('VNTRs/TRF/{sample}.{{chromosome}}.TR.bed',sample=get_reference_ID()),
        beds = expand('graphs/minigraph/{{chromosome}}.{sample}.bed',sample=pangenome_samples)
    output:
        'VNTRs/minigraph/{chromosome}.overlaps.bed'
    shell:
        '''
        gfatools bubble {input.gfa} | cut -f -3,12 |\
        bedtools intersect -f 0.0000000001 -wo -b - -a {input.TRs} |\
        bedtools merge -c 4,8 -o distinct,collapse |\
        python {workflow.basedir}/mg.py {input.beds} > {output}
        '''

rule convert_minigraph:
    input:
        TRs = 'VNTRs/minigraph/{chromosome}.overlaps.bed',
        beds = expand('graphs/minigraph/{{chromosome}}.{sample}.bed',sample=pangenome_samples),
    output:
        'VNTRs/minigraph/{chromosome}.liftover.bed'
    run:
        with open(input.TRs) as fin:
            for line in fin:
                *TR_information, bubbles = line.split()
                source, *_, sink = bubbles.split(',')
                print('\t'.join(TR_information),end='\t',file=output[0])

                TR_coordinates = dict()

                for node,key in zip((source,sink),(4,5)):
                    for hit in subprocess.run(f"awk '${key}==\">{node}\"' {input.beds}",shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]:
                        hit = hit.split()[-1]
                        if hit == '.':
                            continue
                        parts = hit.split(':')
                        if parts[3] in TR_coordinates:
                            TR_coordinates[parts[3]] += f'-{parts[key]},{parts[2]}'
                        elif key == 4:
                            TR_coordinates[parts[3]] = parts[key]

                print('\t'.join((f'{sample}:{coord}' for sample,coord in TR_coordinates.items())),file=output[0])



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
