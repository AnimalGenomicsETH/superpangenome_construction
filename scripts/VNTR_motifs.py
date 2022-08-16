#!/usr/bin/env python
breeds = ['Angus', 'Bison', 'Brahman', 'BSW', 'Gaur', 'Highland', 'Nellore', 'OBV', 'Pied', 'Simmental', 'UCD', 'Yak']

inverted_chrs = {'Angus':['1', '2', '4', '7', '9', '13', '14', '18', '20', '23', '24', '29'], 'Brahman': ['4', '5', '7', '13', '15', '17', '18', '22', '25', '26', '29']}

config = {}

import regex
import subprocess
from collections import defaultdict
def count_VNTRs(sequences,TR):
    allowed_errors = int(len(TR)*0.25)
    counts = dict()
    TR_motifs = dict()
    for asm, sequence in sequences.items():
        asm_ID = asm.split('_')[1].split(':')[0] if '_' in asm else asm
        try:
            temp_motifs = regex.findall(f"(?eV1)({TR}){{2i+2d+1s<={allowed_errors}}}",sequence,concurrent=True,timeout=config.get('timeout',60))
            counts[asm_ID] = len(temp_motifs)
            TR_motifs[asm_ID] = temp_motifs
        except TimeoutError:
            counts[asm_ID] = -1
            break
    return counts, TR_motifs

def get_flank_size(TR_length):
    return 10*TR_length

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
        ix,low,high = name.split(':')[3:6]
        ch,asm = ix[:ix.index('_')],ix[ix.index('_')+1:]
        low, high = int(low), int(high)

        if (ch,asm) in region_per_asm:
            if low < region_per_asm[(ch,asm)][0]:
                region_per_asm[(ch,asm)][0] = low
            if high > region_per_asm[(ch,asm)][1]:
                region_per_asm[(ch,asm)][1] = high
        else:
            region_per_asm[(ch,asm)] = [low,high]

    for (ch,asm),(low,high) in region_per_asm.items():
        offset = get_flank_size(TR_length)
        header, sequence = subprocess.run(f'samtools faidx /cluster/work/pausch/danang/psd/scratch/assembly/{ch}/{asm}_{ch}.fa {ch}_{asm}:{low-offset}-{high+offset} | seqtk seq -l 0 -U {"-r" if ch in inverted_chrs.get(asm,[]) else ""}',shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]
        sequences[header] = sequence
    return sequences

from collections import Counter
def process_line_M(line):
    chrom, start, end, TR_start, TR_end, nodes, TR = line.rstrip().split()
    count = nodes.count(',')+1
    if True:

        regions = []
        for node in nodes.split('|'):
            source, *_, sink = node.split(',')
            regions.extend(subprocess.run(f"awk '$4==\">{source}\"&&$5==\">{sink}\"' {chrom}/*bed",shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1])
        sequences = extract_fasta_regions_M(regions,len(TR))
        if len(sequences)/len(breeds) >= config.get('missing_rate',0):    
            #TEMP LINE
            all_motifs = []
            counts, motifs = count_VNTRs(sequences,TR)
            for a,v in motifs.items():
                print(a,Counter(v).most_common(int(sys.argv[2])))
                all_motifs.extend(v)
            print("All",Counter(all_motifs).most_common(int(sys.argv[2])))
    return

import sys
with open(sys.argv[1]) as fin:
    for line in fin:
        process_line_M(line)
