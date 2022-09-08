intersect_param = {
        "25bp":"max_dist=25",
        "lenient":"max_dist_linear=0.1 max_dist=1000",
        "C":"max_dist_linear=1",
        "OM":"max_dist_linear=1 max_dist=1000"
        }

rule all:
    input:
        expand('merged/{chrom}_OM_merged.vcf',chrom=range(1,30))

rule jasmine:
    input:
        expand('{source}/{{chrom}}.vcf',source=['cactus','pggb','minigraph','assembly','optical'])
    output:
        'merged/{chrom}_{dist}_merged.vcf'
    params:
        _input = lambda wildcards, input: ','.join(input),
        reference = '/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa',
        dist = lambda wildcards: intersect_param[wildcards.dist]
    conda:
        'jasmine'
    threads: 4
    resources:
        mem_mb= 5000,
        walltime= '4:00',
        disk_scratch = 5,
    shell:
        '''
        jasmine --comma_filelist file_list={params._input} threads={threads} out_file={output} out_dir=$TMPDIR \
        genome_file={params.reference} --pre_normalize --dup_to_ins --ignore_strand --allow_intrasample --normalize_type --keep_var_ids \
        {params.dist}
        '''
