include: 'utility.py'

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
        bcftools view --threads {threads} -s {params.samples} -c 1 -i '(ALT="<DEL>"||ALT="<INS>"||ALT="<DUP>")&&ABS(INFO/SVLEN)<1000000' {output.all} | bcftools norm -d none | bcftools annotate -x ^INFO/SVLEN,^INFO/SVTYPE,^FORMAT/GT | bcftools +scatter -o {params._dir} -s {params.chromosomes}
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
        'vcfs/jasmine/{chromosome}.{_group,calls|all}.{setting,lenient|optical|strict|strictest|lenientest|stricter}.vcf'
    params:
        _input = lambda wildcards, input: ','.join(input.vcfs),
        settings = lambda wildcards: config['intersection_parameters'][wildcards.setting]
    conda:
        'jasmine'
    threads: 1
    resources:
        mem_mb= 10000,
        walltime= '4:00',
        scratch = '5G'
    shell:
        '''
        java -Xmx6048m -jar /cluster/work/pausch/alex/software/Jasmine/jasmine.jar \
        --comma_filelist file_list={params._input} threads={threads} out_file={output} out_dir=$TMPDIR \
        genome_file={input.reference} --pre_normalize --ignore_strand --allow_intrasample --normalize_type \
        {params.settings}
        '''

rule stratify_jasmine_repeats:
    input:
        vcf = rules.jasmine.output,
        repeats = expand('assemblies/{chromosome}/{sample}.fa_rm.bed',sample=get_reference_ID(),allow_missing=True)
    output:
        non_repetitive = 'vcfs/jasmine/{chromosome}.{_group,calls|all}.{setting}.non_repetitive.vcf',
        repetitive = 'vcfs/jasmine/{chromosome}.{_group}.{setting}.repetitive.vcf'
    resources:
        walltime = '10'
    params:
        overlap = 0.5
    shell:
        '''
        grep -v "SUPP_VEC=00001" {input.vcf} |\
        bedtools subtract -A -f {params.overlap} -a - -b {input.repeats} > {output.non_repetitive}
        grep -v "SUPP_VEC=00001" {input.vcf} |\
        bedtools intersect -u -f {params.overlap} -a - -b {input.repeats} > {output.repetitive}
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
        mawk 'length($3)==1&&length($4)==1 {{SNP[$5]+=1;next}} {{ if ($3~/,/||$4~/,/) {{MULTI[$5]+=1}} else {{INDEL[$5]+=1}} }} END {{for (key in SNP) {{ print "SNP",key,SNP[key]}} for (key in INDEL) {{ print "INDEL",key,INDEL[key] }} for (key in MULTI) {{ print "MULTI",key,MULTI[key]}} }}' {input} > {output}
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
        echo "pangenome variant sample count" > {output}
        cat {input} >> {output}
        '''

localrules: AF_breakdown
rule AF_breakdown:
    input:
        SVs = expand('vcfs/{{pangenome}}/{chromosome}.SV.vcf',chromosome=range(1,30))
    output:
        'vcfs/{pangenome}/AFs.csv'
    shell:
        '''
        bcftools concat {input} | bcftools query -f '%AC\n' | mawk '{{ A[$1]+=1 }} END {{ for (key in A) {{ print "{wildcards.pangenome}",key,A[key] }} }}' > {output}
        '''

localrules: gather_AF
rule gather_AF:
    input:
        expand('vcfs/{pangenome}/AFs.csv',pangenome=('minigraph','pggb','cactus','assembly'))
    output:
        'vcfs/AFs.csv'
    shell:
        '''
        echo "pangenome AF count" > {output}
        cat {input} >> {output}
        '''

localrules: SV_breakdown
rule SV_breakdown:
    input:
        expand('vcfs/{pangenome}/{chromosome}.SV.vcf',chromosome=range(1,30),pangenome=('minigraph','cactus','pggb','assembly'))
    output:
        'vcfs/SV_breakdown.csv'
    shell:
        '''
        echo "pangenome chromosome deletions insertions other" > {output}
        for p in pggb cactus minigraph assembly
        do
          for c in {{1..29}}
          do
            echo -n "$p $c "
            echo $(bcftools stats vcfs/$p/$c.SV.vcf |\
            awk ' {{ if ($5~/others/) {{ O+=$6 }} else {{ if ($1=="IDD") {{ if ($3~/-/) {{ D+=$4 }} else {{ I+=$4 }} }} }} }} END {{ print D,I,O}} ')
          done
        done >> {output}
'''

rule path_conflicts:
    input:
        expand('vcfs/{pangenome}/{chromosome}.raw.vcf',chromosome=range(1,30),allow_missing=True)
    output:
        'vcfs/{pangenome}/conflicts.txt'
    params:
        samples = list(pangenome_samples.keys()),
        vcf = '$TMPDIR/concat.vcf'
    shell:
        '''
        if [ {wildcards.pangenome} = "minigraph" ]
        then
          for i in {params.samples}; do echo "minigraph SV" $i  $(grep -h "\." graphs/minigraph/*.$i.bed | wc -l) $(cat graphs/minigraph/*.$i.bed | wc -l); done > {output}
        else
          bcftools concat -o $TMPDIR/concat.vcf {input}
          bcftools view -H -i 'abs(ILEN)>=50' -o $TMPDIR/SV.vcf {params.vcf}
          grep -hoP "CONFLICT=\K[A-Z,]+" $TMPDIR/SV.vcf | tr ',' '\\n' | mawk -v T=$(wc -l $TMPDIR/SV.vcf | cut -d' ' -f 1) '{{ A[$1]+=1 }} END {{ for (key in A) {{ print "{wildcards.pangenome}","SV",key,A[key],T }} }}' > {output}
          bcftools view -H -e 'abs(ILEN)>=50' -o $TMPDIR/small.vcf {params.vcf} 
          grep -hoP "CONFLICT=\K[A-Z,]+" $TMPDIR/small.vcf | tr ',' '\\n' | mawk -v T=$(wc -l $TMPDIR/small.vcf | cut -d' ' -f 1) '{{ A[$1]+=1 }} END {{ for (key in A) {{ print "{wildcards.pangenome}","small",key,A[key],T }} }}' >> {output}
        fi
        '''

localrules: gather_conflicts
rule gather_conflicts:
    input:
        expand('vcfs/{pangenomes}/conflicts.txt',pangenomes=('pggb','cactus','minigraph'))
    output:
        'vcfs/conflicts.txt'
    shell:
        '''
        echo "pangenome variant sample count total" > {output}
        cat {input} >> {output}
        '''
