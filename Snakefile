from pathlib import Path,PurePath




pangenome_samples = config['pangenome_samples']
additional_samples = config['additional_samples']
wildcard_constraints:
    chromosome = r'\d+',
    pangenome = r'pggb|cactus|minigraph|assembly'

include: 'snakepit/utility.py'
include: 'snakepit/construct_pangenomes.smk'
include: 'snakepit/decompose_pangenomes.smk'
include: 'snakepit/edit_distance.smk'
include: 'snakepit/VNTRs.smk'
include: 'snakepit/vcf_analysis.smk'

rule all:
    input:
        #graph pangenomes
        expand('graphs/{pangenome}/{chromosome}.gfa',pangenome=('minigraph','pggb','cactus'),chromosome=range(1,30)),
        #variants and their overlaps
        expand('vcfs/{pangenome}/{chromosome}.SV.vcf',pangenome=('minigraph','pggb','cactus','assembly'),chromosome=range(1,30)),
        expand('vcfs/jasmine/{chromosome}.{_group}.{setting}.vcf',chromosome=range(1,30),_group='calls',setting=('strict','lenient')),
        expand('vcfs/jasmine/{chromosome}.{_group}.{setting}.vcf',chromosome=range(1,30),_group='all',setting=('optical',)),
        #edit distances
        expand('edit_distance/{sample}.{chromosome}.{pangenome}.{preset}.{trimmed}.stat',sample=pangenome_samples,chromosome=range(1,30),pangenome=('minigraph','pggb','cactus'),trimmed=('trimmed','untrimmed'),preset='lenient'),
        expand('edit_distance/{sample}.{chromosome}.{pangenome}.{preset}.{trimmed}.stat',sample=additional_samples,chromosome=range(1,30),pangenome=('minigraph','pggb','cactus'),trimmed=('trimmed','untrimmed'),preset='lenient')
