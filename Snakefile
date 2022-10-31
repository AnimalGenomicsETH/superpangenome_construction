pangenome_samples = config['pangenome_samples']
additional_samples = config['additional_samples']
wildcard_constraints:
    chromosome = r'\d+',
    pangenome = r'pggb|cactus|minigraph|assembly'

include: 'snakepit/construct_pangenomes.smk'
include: 'snakepit/decompose_pangenomes.smk'
include: 'snakepit/edit_distance.smk'
#include: 'snakepit/VNTRs.smk'

def get_reference_ID():
    for sample, preset in config['pangenome_samples'].items():
        if preset == 'reference':
            return sample

reference_ID = get_reference_ID()

rule all:
    input:
        expand('graphs/{pangenome}/{chromosome}.gfa',pangenome=('minigraph','pggb','cactus'),chromosome=range(1,30)),
        expand('vcfs/{pangenome}/{chromosome}.SV.vcf',pangenome=('minigraph','pggb','cactus','assembly'),chromosome=range(1,30)),
        #jasmine stuff
        expand('edit_distance/{sample}.{chromosome}.{pangenome}.{trimmed}.stat',sample=pangenome_samples,chromosome=range(1,30),pangenome=('minigraph','pggb','cactus'),trimmed=('trimmed','untrimmed')),
        expand('edit_distance/{sample}.{chromosome}.{pangenome}.{trimmed}.stat',sample=additional_samples,chromosome=range(1,30),pangenome=('minigraph','pggb','cactus'),trimmed=('trimmed','untrimmed'))
