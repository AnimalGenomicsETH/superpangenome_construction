from pathlib import Path,PurePath




pangenome_samples = config['pangenome_samples']
additional_samples = config['additional_samples']
all_samples = list(pangenome_samples.keys()) + list(additional_samples)
wildcard_constraints:
    chromosome = r'\d+',
    pangenome = r'pggb|cactus|minigraph|assembly|miniwaf'

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
        expand('graphs/{pangenome}/stats.yaml',pangenome=('minigraph','pggb','cactus')),
        #variants and their overlaps
        expand('vcfs/{pangenome}/{chromosome}.SV.vcf',pangenome=('minigraph','pggb','cactus','assembly'),chromosome=range(1,30)),
        expand('vcfs/jasmine/{_group}.{setting}.stat',_group='calls',setting=('lenient')),
        expand('vcfs/jasmine/{_group}.{setting}.stat',_group='all',setting=('optical',)),
        expand('vcfs/isec/{mode}.txt',mode=('none','some')),
        #edit distances
        #expand('edit_distance/{sample}.{chromosome}.{pangenome}.{preset}.{trimmed}.gaf',sample=pangenome_samples,chromosome=range(1,30),pangenome=('minigraph','pggb','cactus'),trimmed=('trimmed','untrimmed'),preset='lenient'),
        'edit_distance/summary.tsv',
        #expand('edit_distance/{sample}.{chromosome}.{pangenome}.{preset}.{trimmed}.stat',sample=additional_samples,chromosome=range(1,30),pangenome=('minigraph','pggb','cactus'),trimmed=('trimmed','untrimmed'),preset='lenient'),
        #VNTRs
        expand('VNTRs/{pangenome}/{chromosome}.VNTR.counts',pangenome=('minigraph','pggb','cactus'),chromosome=range(1,30)),
        #Repeats
        expand('graphs/{pangenome}/{chromosome}.fa.masked',pangenome=('minigraph','pggb','cactus'),chromosome=range(1,30))
