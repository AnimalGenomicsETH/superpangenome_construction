breeds = ['Angus', 'Bison', 'Brahman', 'BSW', 'Gaur', 'Highland', 'Nellore', 'OBV', 'Pied', 'Simmental', 'UCD', 'Yak']

#align asm to ref, then determine reversed orders
#sort -k11,11nr angus_reverse.paf | grep -vE "(NKL|X|Y)" | head -n 200| sort -k6,6n | awk '{print $6,$5,$11}' | sort -k2,2h | awk '{seen[$1" "$2]+=$3} END { for (key in seen) { print key,seen[key] } }' |  sort -k1,1n -k3,3nr | awk '!orient[$1]{orient[$1]=$2} END { for (key in orient) { print key,orient[key] } }'
inverted_chrs = {'Angus':['1', '2', '4', '7', '9', '13', '14', '18', '20', '23', '24', '29'], 'Brahman': ['4', '5', '7', '13', '15', '17', '18', '22', '25', '26', '29']}

rule all:
    input:
        'VNTRs.90.csv'

wildcard_constraints:
    chr = r'\d*',
    asm = r'|'.join(breeds),
    rate = r'\d+'

rule minigraph_call:
    input:
        fasta = config['fasta_path'] + '{chr}/{asm}_{chr}.fa',
        gfa = config['gfa_path'] + 'minigraph_{chr}_bs.gfa' 
    output:
        '{chr}/{asm}.bed'
    resources:
        mem_mb = 15000
    shell:
        '''
        minigraph --call -cxasm -t 1 {input.gfa} {input.fasta} > {output}
        '''

rule gfatools_bubble:
    input:
        expand(config['gfa_path'] + 'minigraph_{N}_bs.gfa',N=range(1,30))
    output:
        'minigraph.SVs.bed'
    shell:
        '''
        gfas=({input})
        for i in {{0..28}}; do gfatools bubble ${{gfas[$i]}} | cut -f -3,12 | sed 's/_UCD//g'; done > {output}
        '''


#bedtools merge -d 100 -c 4,5,6 -o distinct,sum,distinct -delim '|'  -i

rule bedtools_intersect:
    input:
        bubbles = 'minigraph.SVs.bed',
        TRs = config['ARS_TR']
    output:
        'putative_VNTRs.{rate}.bed'
    params:
        threshold = 0.000000001 #lambda wildcards: int(wildcards.rate)/100
    shell:
        '''
        bedtools intersect -f {params.threshold} -wo -a {input.bubbles} -b {input.TRs} |\
        awk '{{n=split($4, a, ","); print $1"\\t"$2"\\t"$3"\\t"$6"\\t"$7"\\t"$4"\\t"n"\\t"$8}}' > {output}
        '''

#bedtools intersect -wo -F .01 -a TR/d0_m100.L100.bed -b minigraph.SVs.bed | bedtools merge -i - -c 4,6,7,8 -o distinct,min,max,collapse -delim '|' | awk 'BEGIN{OFS="\t";} {print $1,$5,$6,$2,$3,$7,$4}' > putative_VNTRs.1.merged.bed

rule bedtools_merge:
    input:
        'putative_VNTRs.{rate}.bed'
    output:
        'putative_VNTRs.{rate}.merged.bed'
    params:
        dist = config.get('merging_distance',0)
    shell:
        '''
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

def get_flank_size(TR_length):
    if 'flanking_region' in config:
        offset = config['flanking_region']
    elif 'flank_factor' in config:
        offset = config['flank_factor']*TR_length
    else:
        offset = 0
    return offset

def extract_fasta_regions_M(regions,TR_length):
    sequences = {}
    region_per_asm = {}
    for region in regions:
        name = region.rstrip().split()[5]
        if name == '.':
            continue
        orientation, ix, low, high = name.split(':')[2:6]
        ch,asm = ix[:ix.index('_')],ix[ix.index('_')+1:]
        low, high = int(low), int(high)

        if (ch,asm) in region_per_asm:
            if low < region_per_asm[(ch,asm)][0]:
                region_per_asm[(ch,asm)][0] = low
            if high > region_per_asm[(ch,asm)][1]:
                region_per_asm[(ch,asm)][1] = high
        else:
            region_per_asm[(ch,asm)] = [low,high,orientation]

    for (ch,asm),(low,high,orientation) in region_per_asm.items():
        offset = get_flank_size(TR_length)
        header, sequence = subprocess.run(f'samtools faidx {config["fasta_path"]}{ch}/{asm}_{ch}.fa {ch}_{asm}:{low-offset}-{high+offset} | seqtk seq -l 0 -U {"-r" if orientation=="-" else ""}',shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]
        #header, sequence = subprocess.run(f'samtools faidx {config["fasta_path"]}{ch}/{asm}_{ch}.fa {ch}_{asm}:{low-offset}-{high+offset} | seqtk seq -l 0 -U {"-r" if ch in inverted_chrs.get(asm,[]) else ""}',shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]
        sequences[header] = sequence
    return sequences

def extract_fasta_region_PC(line,TR_length):
    sequence={}
    offset = get_flank_size(TR_length)
    fastadir="/cluster/work/pausch/danang/psd/scratch/assembly"
    for comp in line:
        token=comp.split(",")
        if len(token) > 2:
            chromo,breed,start,stop,strand=comp.split(",")
            if strand == "+":
                start = int(start) - offset
                stop = int(stop) + offset
                seqfa=subprocess.run(f"samtools faidx {fastadir}/{chromo}/{breed}_{chromo}.fa {chromo}_{breed}:{start}-{stop} | seqtk seq -U -l 0",
                        shell=True,capture_output=True).stdout.decode("utf-8").split("\n")[1].upper()
                sequence[breed] = seqfa
            elif strand == "-":
                stop = int(stop) - offset
                start = int(start) + offset
                seqfa=subprocess.run(f"samtools faidx {fastadir}/{chromo}/{breed}_{chromo}.fa {chromo}_{breed}:{stop}-{start} | seqtk seq -U -l 0 -r",
                                    shell=True,capture_output=True).stdout.decode("utf-8").split("\n")[1].upper()
                sequence[breed] = seqfa
    return sequence


import math
from multiprocessing import Pool

import Levenshtein
#statistics.variance throws an error if only 1 sample, so use numpy
#from statistics import variance
from numpy import var as variance
def postprocess_line(chrom,start,end,TR_start,TR_end,TR,sequences,counts):
    stats = []
    for B in breeds:
        stats.append(counts.get(B,math.nan))
    return ','.join(map(str,[chrom,start,end,TR_start,TR_end,TR,len(sequences)] + stats))


def process_line_M(line):
    chrom, start, end, TR_start, TR_end, nodes, TR = line.rstrip().split()
    count = nodes.count(',')+1
    if '|' in TR:
        TR = Levenshtein.median(TR.split('|'))
    if config.get('low',2) < count < config.get('high',math.inf) and config.get('TR_low',5) <= len(TR) <= config.get('TR_high',math.inf):

        regions = []
        for node in nodes.split('|'):
            source, *_, sink = node.split(',')
            regions.extend(subprocess.run(f"awk '$4==\">{source}\"&&$5==\">{sink}\"' {chrom}/*bed",shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1])
        sequences = extract_fasta_regions_M(regions,len(TR))
        if len(sequences)/len(breeds) >= config.get('missing_rate',0):    
            #TEMP LINE
            counts = count_VNTRs(sequences,TR)
            return postprocess_line(chrom,start,end,TR_start,TR_end,TR,sequences,counts)
    return

def process_line_PC(line):
    try:
        chrom, start, end, TR, *regions = line.rstrip().split()
        if config.get('TR_low',5) <= len(TR) <= config.get('TR_high',math.inf):
            sequences = extract_fasta_region_PC(regions,len(TR))

            counts = count_VNTRs(sequences,TR)
            return postprocess_line(chrom,start,end,start,end,TR,sequences,counts)
    except:
        pass
    return

from itertools import product
rule process_M_VNTRs:
    input:
        VNTRs = 'putative_VNTRs.{rate}.merged.bed',
        beds = expand('{chr}/{asm}.bed',chr=range(1,30),asm=breeds),
        fasta = expand(config['fasta_path'] + '{chr}/{asm}_{chr}.fa',chr=range(1,30),asm=breeds)
    output:
        'VNTRs.{rate}.csv'
    threads: 32
    resources:
        mem_mb = 200,
        walltime = '24:00'
    run:
        with open(input.VNTRs,'r') as fin, open(output[0],'w') as fout:
            print('chr,start,end,TR_start,TR_end,TR,alignments,' + ','.join(breeds),file=fout)
            
            if config.get('debug',False):
                for i,line in enumerate(fin):
                    result = process_line_M(line)
                    if result:
                        print(result,file=fout)
                        fout.flush()
            else:
                with Pool(threads) as p:
                    for i, result in enumerate(p.imap_unordered(process_line_M, fin, 50)):
                        if i % 100 == 0:
                            print(f'done {i=} results',flush=True)
                            fout.flush()
                        if result:
                            print(result,file=fout)

rule process_regions_PC:
    input:
        regions = 'TR_tsvs/{pangenome}.{chr}.tsv'
    output:
        'TR_tsvs/{pangenome}.{chr}.count'
    threads: 8
    resources:
        mem_mb = 500,
        walltime = "24:00"
    run:
        with open(input.regions,'r') as fin, open(output[0],'w') as fout:
            print('chr,start,end,TR_start,TR_end,TR,alignments,' + ','.join(breeds),file=fout)
            if config.get('debug',False):
                pass
            else:
                with Pool(threads) as p:
                    for i, result in enumerate(p.imap_unordered(process_line_PC, fin, 50)):
                        if i % 100 == 0:
                            print(f'done {i=} results',flush=True)
                            fout.flush()
                        if result:
                            print(result,file=fout)
