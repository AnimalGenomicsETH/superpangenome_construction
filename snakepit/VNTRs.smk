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

rule gfatools_bubble:
    input:
        gfa = 'graphs/minigraph/{chromosome}.basic.gfa',
        TRs = expand('VNTRs/TRF/{sample}.{{chromosome}}.TR.bed',sample=get_reference_ID()),
        beds = expand('graphs/minigraph/{{chromosome}}.{sample}.bed',sample=pangenome_samples)
    output:
        'VNTRs/minigraph/{chromosome}.liftover.bed'
    shell:
        '''
        gfatools bubble {input.gfa} | cut -f -3,12 |\
        bedtools intersect -f 0.0000000001 -wo -b - -a {input.TRs} |\
        bedtools merge -c 4,8 -o distinct,collapse |\
        python {workflow.basedir}/minigraph_VNTR_conversion.py {input.beds} > {output}
        '''

import regex
punctuation = regex.compile(r'[:,-]')

def extract_fasta(regions,chromosome,TR_length):
    sequences = {}
    offset = 500#get_flank_size(TR_length)
    for region in regions:
        sample, start, stop, orientation = punctuation.split(region)
        if orientation == '+':
            start = int(start) - offset
            stop = int(stop) + offset
            extra_flag = ''
        else:
            start = int(stop) - offset
            stop = int(start) + offset
            extra_flag = '-r'
        sequences[sample] = subprocess.run(f'samtools faidx assemblies/{chromosome}/{sample}.fa {sample}:{start}-{stop} | seqtk seq -U -l 0 {extra_flag}',shell=True,capture_output=True).stdout.decode("utf-8").split("\n")[1]

    return sequences

import math
from collections import defaultdict
def count_VNTRs(sequences,TR):
    allowed_errors = int(len(TR)*config.get('TR_divergence_limit',0))# if config.get('scaling','linear') != 'log' else math.floor(math.sqrt(len(TR)))
    counts = dict()
    for sample, sequence in sequences.items():
        try:
            counts[sample] = len(regex.findall(f"(?eV1)({TR}){{2i+2d+1s<={allowed_errors}}}",sequence,concurrent=True,timeout=config.get('timeout',60)))
        except TimeoutError:
            counts[asm_ID] = -1
            break
    return counts

def process_VNTR_line(line,chromosome,samples):
    if True:
        chrom, start, end, TR, *regions = line.rstrip().split()
        if config.get('TR_low',5) <= len(TR) <= config.get('TR_high',math.inf):
            sequences = extract_fasta(regions,chromosome,len(TR))

            counts = count_VNTRs(sequences,TR)
            return ','.join(map(str,(chromosome,start,end,TR,len(sequences),*[counts.get(sample,math.nan) for sample in samples])))
            return postprocess_line(chrom,start,end,start,end,TR,sequences,counts)
    #except Exception as e:
    #    print(e)
    return

from functools import partial

rule process_VNTRs:
    input:
        'VNTRs/{pangenome}/{chromosome}.liftover.bed'
    output:
        'VNTRs/{pangenome}/{chromosome}.VNTR.counts'
    params:
        samples = pangenome_samples
    run:
        with open(input[0],'r') as fin, open(output[0],'w') as fout:
            #print('chr,start,end,TR_start,TR_end,TR,alignments,' + ','.join(),file=fout)
            
            if config.get('debug',False):
                for i,line in enumerate(fin):
                    result = process_VNTR_line(line,wildcards.chromosome,params.samples)
                    if result:
                        print(result,file=fout)
                        fout.flush()
            else:
                with Pool(threads) as p:
                    for i, result in enumerate(p.imap_unordered(partial(process_VNTR_line,chromosome=wildcards.chromosome,samples=params.samples), fin, 50)):
                        if i % 100 == 0:
                            print(f'done {i=} results',flush=True)
                            fout.flush()
                        if result:
                            print(result,file=fout)


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
