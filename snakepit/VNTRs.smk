

rule odgi_position:
    input:
        gfa = 'graphs/{pangenome}/{chromosome}.gfa',
        TRs = 'VNTRs/TRF/{chromosome}.TRs.bed'
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
        gfas=({input})
        for i in {{0..28}}; do gfatools bubble ${{gfas[$i]}} | cut -f -3,12 | sed 's/_UCD//g'; done > {output}
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
