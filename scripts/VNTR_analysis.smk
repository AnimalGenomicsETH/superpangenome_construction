
breeds = ['Angus', 'Bison', 'Brahman', 'BSW', 'Gaur', 'Highland', 'Nellore', 'OBV', 'Pied', 'Simmental', 'UCD', 'Yak']
rule all:
    input:
        'VNTRs.90.csv',
        expand('{chr}/{asm}.bed',chr=range(1,30),asm=breeds)

wildcard_constraints:
    chr = r'\d*',
    asm = r'|'.join(breeds)

rule minigraph_call:
    input:
        fasta = '/cluster/work/pausch/danang/psd/scratch/assembly/{chr}/{asm}_{chr}.fa',
        gfa = '/cluster/work/pausch/to_share/sv_analysis/new_minigraph/minigraph_{chr}_bs.gfa' 
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
        #expand(config['gfa_path'] + 'minigraph_{N}_bs.gfa',N=range(1,30))
        expand('/cluster/work/pausch/to_share/sv_analysis/new_minigraph/minigraph_{N}_bs.gfa',N=range(1,30))
    output:
        'minigraph.SVs.bed'
    shell:
        '''
        gfas=({input})
        for i in {{0..28}}; do gfatools bubble ${{gfas[$i]}} | cut -f -3,12 | sed 's/_UCD//g'; done > {output}
        '''

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
        bedtools intersect -f {params.threshold} -wo -a {input.bubbles} -b {input.TRs} | awk '{{n=split($4, a, ","); print $1"\\t"$2"\\t"$3"\\t"$4"\\t"n"\\t"$8}}' > {output}
        '''

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

import regex
def count_VNTRs(sequences,TR):
    allowed_errors = int(len(TR)*config.get('TR_divergence_limit',0))
    return {asm: len(regex.findall(f"({TR}){{e<={allowed_errors}}}",sequence)) for asm, sequence in sequences.items()}

def extract_fasta_regions(regions):
    sequences = {}
    for region in regions:
        parts= region.rstrip().split()
        name = parts[5]
        if name == '.':
            continue
        ix,low,high = name.split(':')[3:6]
        offset = config.get('offset',500)
        ch,asm = ix[:ix.index('_')],ix[ix.index('_')+1:]
        header, sequence = subprocess.run(f'samtools faidx /cluster/work/pausch/danang/psd/scratch/assembly/{ch}/{asm}_{ch}.fa {ix}:{int(low)-offset}-{int(high)+offset} | seqtk seq -l 0',shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]
        
        sequences[header] = sequence.upper() if 'Angus' not in header else reverse_complement(sequence.upper())
    return sequences

import math
rule process_VNTRs:
    input:
        VNTRs = 'putative_VNTRs.{rate}.bed',
        beds = expand('{chr}/{asm}.bed',chr=range(1,30),asm=breeds),
        fasta = expand('/cluster/work/pausch/danang/psd/scratch/assembly/{chr}/{asm}_{chr}.fa',chr=range(1,30),asm=breeds)
    output:
        'VNTRs.{rate}.csv'
    run:
        with open(input.VNTRs,'r') as fin, open(output[0],'w') as fout:
            print('chr,start,end,TR,' + ','.join(breeds),file=fout)
            for line in fin:
                chrom, start, end, nodes, count, TR = line.rstrip().split()
                count = int(count)
                if config.get('low',2) < count < config.get('high',math.inf) and config.get('TR_low',5) <= len(TR) <= config.get('TR_high'):
                    source, *_, sink = nodes.split(',')
                    regions = subprocess.run(f"awk '$4==\">{source}\"&&$5==\">{sink}\"' {chrom}/*bed",shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]
                    sequences = extract_fasta_regions(regions)

                    counts = count_VNTRs(sequences,TR)
                    if len(counts)/len(breeds) < config.get('missing_rate',0):
                        continue
                    counts_clean = {regex.match( r'.*?_(.*):.*',k).group(1):v for k,v in counts.items()}
                    print(','.join(map(str,[chrom,start,end,TR] + [counts_clean.get(b,math.nan) for b in breeds])),file=fout)
                    
        #process lines, select good ones
        #extract local fasta
        #count VNTRs

#rule localise_bed:
#    input:
#        expand('{{chr}}/{asm}.bed',asm=breeds)
#    output:
#        ''
#    shell:
#        '''
#        awk '$1==S&&$2==$3'
#        '''
