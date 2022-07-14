
rule all:
    input:
        'TR/ARS.TR.nonrep.bed'

rule TRF:
    output:
        'TR/{asm}.fa.2.7.7.80.10.50.500.dat'
        #'input/{asm}.TR.bed'
    threads: 1
    resources:
        mem_mb = 1000,
        walltime = '4:00'
    params:
        lambda wildcards, output: ' '.join(output[0].split('.')[-8:-1])
    shell:
        '''
        TRF {config[reference]} {params} -h
        '''

rule postprocess_TRF:
    input:
        'TR/{asm}.fa.2.7.7.80.10.50.500.dat'
        #'input/{asm}.TR.bed'
    output:
        'TR/{asm}.TR.filtered.bed'
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
        python postprocess_TRs.py > {output}
        '''


rule finalise_TR:
    input:
        TRs = 'TR/{asm}.TR.filtered.bed',
        repeats = config['repeats']
    output:
        ctrl = 'TR/{asm}.TR.ctrl.bed',
        TR = 'TR/{asm}.TR.nonrep.bed',
    shell:
        '''
        #no overlaps -> keep "good" TRs
        bedtools subtract -A -f 0.9 -F 0.000001 -a {input.TRs} -b {input.repeats} > {output.TR}
        #finds overlaps -> repetitive TRs so filter
        bedtools intersect -u -f 0.9 -F 0.000001 -a {input.TRs} -b {input.repeats} > {output.ctrl}
        '''
