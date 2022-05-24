breeds = ['Angus', 'Bison', 'Brahman', 'BSW', 'Gaur', 'Highland', 'Nellore', 'OBV', 'Pied', 'Simmental', 'UCD', 'Yak']

#sort -k11,11nr angus_reverse.paf | grep -vE "(NKL|X|Y)" | head -n 200| sort -k6,6n | awk '{print $6,$5,$11}' | sort -k2,2h | awk '{seen[$1" "$2]+=$3} END { for (key in seen) { print key,seen[key] } }' |  sort -k1,1n -k3,3nr | awk '!orient[$1]{orient[$1]=$2} END { for (key in orient) { print key,orient[key] } }'
inverted_chrs = {'Angus':['1', '2', '4', '7', '9', '13', '14', '18', '20', '23', '24', '29'], 'Brahman': ['4', '5', '7', '13', '15', '17', '18', '22', '25', '26', '29']}

rule all:
    input:
        'VNTRs.90.csv'

wildcard_constraints:
    chr = r'\d*',
    asm = r'|'.join(breeds)

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
        threshold = lambda wildcards: int(wildcards.rate)/100
    shell:
        '''
        #filters out duplicate lines where TRF can output multiple TRs
        bedtools intersect -f {params.threshold} -wo -a {input.bubbles} -b {input.TRs} | awk '{{n=split($4, a, ","); print $1"\\t"$2"\\t"$3"\\t"$4"\\t"n"\\t"$8}}' | awk '!seen[$1$2$3$4]++' > {output}
        '''

import regex
from collections import defaultdict
def count_VNTRs(sequences,TR):
    allowed_errors = int(len(TR)*config.get('TR_divergence_limit',0))
    #could implement some sqrt approach, but then unequal meaning of divergence
    counts = defaultdict(int)
    for asm, sequence in sequences.items():
        counts[asm.split('_')[1].split(':')[0]] += len(regex.findall(f"({TR}){{e<={allowed_errors}}}",sequence,concurrent=True))

    return dict(counts)
    #return {asm.split('_')[1].split(':')[0]: len(regex.findall(f"({TR}){{e<={allowed_errors}}}",sequence,concurrent=True)) for asm, sequence in sequences.items()}

def extract_fasta_regions(regions):
    sequences = {}
    for region in regions:
        parts= region.rstrip().split()
        name = parts[5]
        if name == '.':
            continue
        ix,low,high = name.split(':')[3:6]
        offset = config.get('flanks',100)
        ch,asm = ix[:ix.index('_')],ix[ix.index('_')+1:]
        
        if (int(high) - int(low)) > config.get('max_region_length',10000):
            print('TOO LONG SKIPPING',region)
            continue

        header, sequence = subprocess.run(f'samtools faidx {config["fasta_path"]}{ch}/{asm}_{ch}.fa {ix}:{int(low)-offset}-{int(high)+offset} | seqtk seq -l 0 -U {"-r" if ch in inverted_chrs.get(asm,[]) else ""}',shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]
        sequences[header] = sequence
    return sequences

import math
from multiprocessing import Pool

import Levenshtein

def process_line(line):
    chrom, start, end, nodes, count, TR = line.rstrip().split()
    count = int(count)
    if config.get('low',2) < count < config.get('high',math.inf) and config.get('TR_low',5) <= len(TR) <= config.get('TR_high',math.inf):
        if '|' in TR:
            TR = Levenshtein.median(TR.split('|'))

        regions = []
        for node in nodes.split('|'):
            source, *_, sink = node.split(',')
            #someway to handle multi-bubble
            regions.extend(subprocess.run(f"awk '$4==\">{source}\"&&$5==\">{sink}\"' {chrom}/*bed",shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1])
        sequences = extract_fasta_regions(regions)

        if len(sequences)/len(breeds) >= config.get('missing_rate',0):
            counts = count_VNTRs(sequences,TR)
            return ','.join(map(str,[chrom,start,end,TR] + [counts.get(b,math.nan) for b in breeds]))
    return

rule process_VNTRs:
    input:
        VNTRs = 'putative_VNTRs.{rate}.bed',
        beds = expand('{chr}/{asm}.bed',chr=range(1,30),asm=breeds),
        fasta = expand(config['fasta_path'] + '{chr}/{asm}_{chr}.fa',chr=range(1,30),asm=breeds)
    output:
        'VNTRs.{rate}.csv'
    threads: 18
    resources:
        mem_mb = 500,
        walltime = '120:00'
    run:
        with open(input.VNTRs,'r') as fin, open(output[0],'w') as fout:
            print('chr,start,end,TR,' + ','.join(breeds),file=fout)
            #for i,result in enumerate([process_line(line) for line in fin]):
            #    print(result,file=fout)

            with Pool(threads) as p:
                for i, result in enumerate(p.imap_unordered(process_line, fin, 10)):
                    if i % 100 == 0:
                        print(f'done {i=} results',flush=True)
                        fout.flush()
                    if result:
                        print(result,file=fout)

            #for line in fin:
            #    chrom, start, end, nodes, count, TR = line.rstrip().split()
            #    count = int(count)
            #    if config.get('low',2) < count < config.get('high',math.inf) and config.get('TR_low',5) <= len(TR) <= config.get('TR_high'):
            #        source, *_, sink = nodes.split(',')
            #        regions = subprocess.run(f"awk '$4==\">{source}\"&&$5==\">{sink}\"' {chrom}/*bed",shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]
            #        sequences = extract_fasta_regions(regions)

            #        counts = count_VNTRs(sequences,TR)
            #        if len(counts)/len(breeds) < config.get('missing_rate',0):
            #            continue
            #        counts_clean = {regex.match( r'.*?_(.*):.*',k).group(1):v for k,v in counts.items()}
            #        print(','.join(map(str,[chrom,start,end,TR] + [counts_clean.get(b,math.nan) for b in breeds])),file=fout)
