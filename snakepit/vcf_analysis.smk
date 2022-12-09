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
    params:
        reference_ID = get_reference_ID()
    shell:
        '''
        bcftools annotate --rename-chrs <(echo -e "{wildcards.chromosome}\t{params.reference_ID}") -o {output} {input}
        '''

rule jasmine:
    input:
        vcfs = lambda wildcards: expand('vcfs/{pangenome}/{{chromosome}}.SV.vcf',pangenome=get_variants(wildcards._group)),
        reference = expand('assemblies/{{chromosome}}/{ref_ID}.fa',ref_ID=get_reference_ID())
    output:
        'vcfs/jasmine/{chromosome}.{_group,calls}.{setting,lenient}.vcf'
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

rule stratify_jasmine_repeats:
    input:
        vcf = rules.jasmine.output,
        repeats = expand('assemblies/{chromosome}/{sample}.fa_rm.bed',sample=get_reference_ID(),allow_missing=True)
    output:
        non_repetitive = 'vcfs/jasmine/{chromosome}.{_group,calls}.{setting}.non_repetitive.vcf',
        repetitive = 'vcfs/jasmine/{chromosome}.{_group}.{setting}.repetitive.vcf'
    resources:
        walltime = '10'
    params:
        overlap = 0.75
    shell:
        '''
        bedtools subtract -A -f {params.overlap} -a {input.vcf} -b {input.repeats} > {output.non_repetitive}
        bedtools intersect -u -f {params.overlap} -a {input.vcf} -b {input.repeats} > {output.repetitive}
        '''

localrules: summarise_jasmine
rule summarise_jasmine:
    input:
        non_repetitive = expand('vcfs/jasmine/{chromosome}.{_group}.{setting}.non_repetitive.vcf',chromosome=range(1,30),allow_missing=True),
        repetitive = expand('vcfs/jasmine/{chromosome}.{_group}.{setting}.repetitive.vcf',chromosome=range(1,30),allow_missing=True)
    output:
        'vcfs/jasmine/{_group}.{setting}.stat'
    shell:
        '''
        grep -hoP "SUPP_VEC=\K\d+" {input.non_repetitive} | mawk '{{A[$1]+=1}} END {{ for (k in A) {{ print "non_repetitive",k,A[k] }} }}' > {output}
        grep -hoP "SUPP_VEC=\K\d+" {input.repetitive} | mawk '{{A[$1]+=1}} END {{ for (k in A) {{ print "repetitive",k,A[k] }} }}' >> {output}
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
        expand(rules.bcftools_isec.output,chromosome=range(1,30),allow_missing=True)
    output:
        'vcfs/isec/{mode}.txt'
    shell:
        '''
        mawk 'length($3)==1&&length($4)==1 {{SNP[$5]+=1;next}} {{INDEL[$5]+=1}} END {{for (key in SNP) {{ print "SNP",key,SNP[key]}} for (key in INDEL) {{ print "INDEL",key,INDEL[key] }} }}' {input} > {output}
        '''


rule sample_breakdown_variants:
    input:
        SVs = expand('vcfs/{{pangenome}}/{chromosome}.SV.vcf',chromosome=range(1,30)),
        small = expand('vcfs/{{pangenome}}/{chromosome}.small.vcf.gz',chromosome=range(1,30))
    output:
        'vcfs/{pangenome}/per_sample.csv'
    params:
        tmp_small = '$TMPDIR/small.vcf.gz',
        key = lambda wildcards: ('$13','$13') if wildcards.pangenome != 'assembly' else ('$9','$5')
    shell:
        '''
        bcftools concat {input.SVs} | bcftools stats -s- | awk '$1=="PSC" {{ print "{wildcards.pangenome}","SV",$3,{params.key[0]} }}' > {output}
        bcftools concat --naive-force -o {params.tmp_small} {input.small}
        bcftools view -i 'TYPE="snp"' {params.tmp_small} | bcftools stats -s- | awk '$1=="PSC" {{ print "{wildcards.pangenome}","SNP",$3,{params.key[1]} }}' >> {output}
        bcftools view -i 'TYPE!="snp"' {params.tmp_small} | bcftools stats -s- | awk '$1=="PSC" {{ print "{wildcards.pangenome}","Indel",$3,{params.key[0]} }}' >> {output}
        '''

localrules: gather_sample_breakdown
rule gather_sample_breakdown:
    input:
        expand('vcfs/{pangenome}/per_sample.csv',pangenome=('minigraph','pggb','cactus','assembly'))
    output:
        'vcfs/per_sample.csv'
    shell:
        '''
        cat {input} > {output}
        '''
