#include: 'utility.py'

def get_variants(_group):
    match _group:
        case 'calls':
            return ('minigraph','cactus','pggb','assembly')
        case 'base-level':
            return ('cactus','pggb','assembly')
        case _:
            return ('minigraph','cactus','pggb','assembly','optical')

rule prepare_optical_maps:
    output:
        all = temp('vcfs/optical/all.vcf.gz'), 
        chromosomes = temp(expand('vcfs/optical/{chromosome}.vcf',chromosome=range(1,30)))
    params:
        URL = config['optical_map_URL'],
        _dir = lambda wildcards, output: PurePath(output.all).parent,
        samples = 'Hereford_1,Hereford_2,Nelore_N2,Nelore_N4689',
        chromosomes = ','.join(map(str,range(1,30)))
    envmodules:
        'eth_proxy'
    threads: 1
    resources:
        mem_mb = 1500
    shell:
        '''
        wget -O {output.all} {params.URL}
        bcftools view --threads {threads} -s {params.samples} -i '(ALT="<DEL>"||ALT="<INS>"||ALT="<DUP>")&&ABS(INFO/SVLEN)<1000000' {output.all} | bcftools norm -d none | bcftools +scatter -o {params._dir} -s {params.chromosomes}
        '''

localrules: rename_optical_chromosomes
rule rename_optical_chromosomes:
    input:
        'vcfs/optical/{chromosome}.vcf'
    output:
        'vcfs/optical/{chromosome}.SV.vcf'
    shell:
        '''
        bcftools annotate --rename-chrs <(echo -e "{wildcards.chromosome}\tHER") -o {output} {input}
        '''

rule jasmine:
    input:
        vcfs = lambda wildcards: expand('vcfs/{pangenome}/{{chromosome}}.SV.vcf',pangenome=get_variants(wildcards._group)),
        reference = expand('assemblies/{{chromosome}}/{ref_ID}.fa',ref_ID=get_reference_ID())
    output:
        'vcfs/jasmine/{chromosome}.{_group}.{setting}.vcf'
    params:
        _input = lambda wildcards, input: ','.join(input.vcfs),
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
        java -jar /cluster/work/pausch/alex/software/Jasmine/jasmine.jar \
        --comma_filelist file_list={params._input} threads={threads} out_file={output} out_dir=$TMPDIR \
        genome_file={input.reference} --pre_normalize --dup_to_ins --ignore_strand --allow_intrasample --normalize_type \
        {params.settings}
        '''

rule bcftools_isec:
    input:
        expand('vcfs/{pangenome}/{{chromosome}}.small.vcf.gz',pangenome=get_variants('base-level'))
    output:
        'vcfs/isec/{chromosome}.{mode}.isec'
    threads: 2
    resources:
        mem_mb = 1500
    shell:
        '''
        bcftools isec -c {wildcards.mode} --threads {threads} -n +1 -o {output} {input}
        '''

localrules: count_isec_overlaps
rule count_isec_overlaps:
    input:
        rules.bcftools_isec.output
    output:
        'vcfs/isec/{chromosome}.{mode}.txt'
    shell:
        '''
        awk 'length($3)==1&&length($4)==1 {{SNP[$5]+=1;next}} {{INDEL[$5]+=1}} END {{for (key in SNP) {{ print "SNP",key,SNP[key]}} for (key in INDEL) {{ print "INDEL",key,INDEL[key] }} }}' {input} > {output}
        '''
