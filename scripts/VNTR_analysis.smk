breeds = ['Angus', 'Bison', 'Brahman', 'BSW', 'Gaur', 'Highland', 'Nellore', 'OBV', 'Pied', 'Simmental', 'UCD', 'Yak']
rule all:
    input:
        expand('{chr}/{asm}.bed',chr=range(1,30),asm=breeds)

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

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def load_region_fasta(fname):
    sequences = {}
    seq = None
    with open(fname,'r') as fin:
        for lin in fin:
            #grab header line, index at 1 to remove leading '>'
            if seq:    
                seq = line.strip()[1:]
            else:
                #all Angus sequences seem to be other strand, so use the reverse complement so the matching is correct
                sequences[seq] == line.strip().upper() if 'Angus' not in seq else reverse_complement(line.strip().upper())
                seq = None
    return sequences

import regex
def count_VNTRs(sequences,TR,allowed_errors=3):
    for asm, sequence in sequences.items():
        print(asm, len(regex.findall(f"{TR}{{e<={allowed_errors}}}",sequence)))

rule gfatools_bubble:
    input:
        expand(config['gfa_path'] + 'minigraph_{N}_bs.gfa',N=range(1,30))
    output:
        'minigraph.SVs.bed'
    shell:
        '''
        gfas =({input})
        for i in {1..29}; do gfatools bubble ${{gfas[$i]}} | cut -f -3,12 | sed 's/_UCD//g'; done > {output}
        '''

rule bedtools_intersect:
    input:
        bubbles = 'minigraph.SVs.bed',
        TRs = config['ARS_TR']
    output:
        'putative_VNTRs.strict.bed'
    shell:
        '''
        bedtools intersect -f 0.9 -a {input.bubbles} -b {input.TRs} | awk '{{n=split($4, a, ","); print $0"\t"n}}' > {output}
        '''

#rule localise_bed:
#    input:
#        expand('{{chr}}/{asm}.bed',asm=breeds)
#    output:
#        ''
#    shell:
#        '''
#        awk '$1==S&&$2==$3'
#        '''
