#include: 'utility.py'

def get_variants(_group):
    if _group == 'calls':
        return ('minigraph','cactus','pggb','assembly')
    else:
        return ('minigraph','cactus','pggb','assembly','optical')

rule prepare_optical_maps:
    output:
        all = temp('vcfs/optical/all.vcf.gz'), 
        chromosomes = expand('vcfs/optical/{chromosome}.SV.vcf',chromosome=range(1,30))
    params:
        URL = config['optical_map_URL'],
        _dir = lambda wildcards: PurePath(output.all).parent,
        samples = 'Hereford_1,Hereford_2,Nelore_N2,Nelore_N4689',
        chromosomes = ','.join(map(str,range(1,30)))
    envs:
        'eth_proxy'
    threads: 1
    resources:
        mem_mb = 1500
    shell:
        '''
        wget -O {output.all} {params.URL}
        bcftools view --threads {threads} -s {params.samples} -i '(ALT="<DEL>"||ALT="<INS>"||ALT="<DUP>")&&ABS(INFO/SVLEN)<1000000' {output.all} | bcftools norm -d none | bcftools +scatter -o optical -s {params.chromosomes}
        rename .vcf .SV.vcf *.vcf
        '''

rule jasmine:
    input:
        lambda wildcards: expand('vcfs/{pangenome}/{{chromosome}}.SV.vcf',pangenome=get_variants(wildcards._group))
    output:
        'vcfs/jasmine/{chromosome}.{_group}.{setting}.vcf'
    params:
        _input = lambda wildcards, input: ','.join(input),
        reference = config['pangenome_samples'][get_reference_ID()],
        settings = lambda wildcards: config['intersection_parameters'][wildcards.setting]
    conda:
        'jasmine'
    threads: 4
    resources:
        mem_mb= 5000,
        walltime= '4:00',
        scratch = '5G'
    shell:
        '''
        jasmine --comma_filelist file_list={params._input} threads={threads} out_file={output} out_dir=$TMPDIR \
        genome_file={params.reference} --pre_normalize --dup_to_ins --ignore_strand --allow_intrasample --normalize_type --keep_var_ids \
        {params.dist}
        '''
