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
        walltime= '4h',
        scratch = '5G'
    shell:
        '''
        java -Xmx6048m -jar /cluster/work/pausch/alex/software/Jasmine/jasmine.jar \
        --comma_filelist file_list={params._input} threads={threads} out_file={output} out_dir=$TMPDIR \
        genome_file={input.reference} --pre_normalize --ignore_strand --allow_intrasample --normalize_type \
        {params.settings}

        sed -i '20i ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">\\n##INFO=<ID=STRANDS,Number=1,Type=String,Description=""> {output}
        '''

localrules: generate_genome_annotation
rule generate_genome_annotation:
    input:
        masked = expand('assemblies/{chromosome}/{sample}.fa_rm.bed',sample=get_reference_ID(),allow_missing=True),
        TR = expand('VNTRs/TRF/{sample}.{chromosome}.TRF',sample=get_reference_ID(),allow_missing=True),
        low_map = config['low_mappability'],
        fai = expand('assemblies/{chromosome}/{sample}.fa.fai',sample=get_reference_ID(),allow_missing=True)
    output:
        total = 'vcfs/genome_annotation.{chromosome}.bed',
        hard = 'vcfs/genome_annotation.hard.{chromosome}.bed',
        normal = 'vcfs/genome_annotation.normal.{chromosome}.bed',
    shell:
        '''
        {{ awk -v OFS='\\t' '/Satellite/ {{ print $1,$2,$3,"0" }} ' {input.masked} | bedtools merge -d 10000 -c 4 -o distinct; \
        awk -v OFS='\\t' '!/Satellite/ {{ print $1,$2,$3,"2" }} ' {input.masked}; \
        awk -v OFS='\\t' 'NR>1 {{ print "HER",$1,$2,"1" }} ' {input.TR} | bedtools sort -i - | bedtools merge -d 0 -i - -c 4 -o first | bedtools slop -g {input.fai} -i - -b 0; \
        awk -v OFS='\\t' -v c={wildcards.chromosome} '$1==c {{ print "HER",$2,$3,"3"}}' {input.low_map} | bedtools slop -g {input.fai} -i - -b 10; }} |\
        bedtools sort -i - |\
        bedtools merge -d 0 -c 4 -o min | tee {output.hard} |\
        bedtools complement -L -g {input.fai} -i - | awk -v OFS='\\t' ' {{ print $1,$2,$3,"Normal" }} ' >> {output.normal}
        cat {output.hard} {output.normal} |\
        bedtools sort -i - |\
        sed -e 's/0$/Satellite/' -e 's/2$/Repetitive/' -e 's/1$/Tandem repeat/' -e 's/3$/Low mappability/' > {output.total}
        '''


#awk -v OFS='\\t' '{{ print $1,$2,$3,"2" }} ' {input.TR}; \

rule stratify_jasmine_repeats:
    input:
        vcf = rules.jasmine.output,
        annotation = expand(rules.generate_genome_annotation.output['total'],sample=get_reference_ID(),allow_missing=True)
    output:
        'vcfs/jasmine/{chromosome}.{_group}.{setting}.stat'
    resources:
        walltime = '10m'
    params:
        overlap = 0.5
    shell:
        '''
        bcftools annotate -x ^INFO/SUPP_VEC {input.vcf} |\
        bcftools annotate -a {input.annotation} -c CHROM,FROM,TO,CLASS -H '##INFO=<ID=CLASS,Number=1,Type=String,Description="Genomic region">' |\
        grep -oP "(SUPP_VEC=\K\d+|CLASS=\K[^\\t]+)" |paste -d "\\t"  - - |\
        awk -v OFS='\\t' -F '\\t' ' {{ ++A[$2][$1] }} END {{ for (c in A) {{ for (k in A[c]) {{ print c,k,A[c][k] }} }} }}' > {output}
        '''

localrules: summarise_jasmine
rule summarise_jasmine:
    input:
        repetitive = expand(rules.stratify_jasmine_repeats.output[0],chromosome=range(1,30),allow_missing=True)
    output:
        'vcfs/jasmine/{_group}.{setting}.stat'
    shell:
        '''
        awk -v OFS='\\t' -F '\\t' ' {{ A[$1][$2]+=$3 }} END {{ for (c in A) {{ for (k in A[c]) {{ print c,k,A[c][k] }} }} }}' {input} > {output}
        '''

rule bcftools_isec:
    input:
        vcfs = expand('vcfs/{pangenome}/{{chromosome}}.small.vcf.gz',pangenome=get_variants('base-level')),
        annotation = expand(rules.generate_genome_annotation.output['total'],sample=get_reference_ID(),allow_missing=True)
    output:
        'vcfs/isec/{chromosome}.{mode}.isec'
    threads: 2
    resources:
        mem_mb = 1500
    shell:
        '''
        bcftools isec -c {wildcards.mode} --threads {threads} -n +1 {input.vcfs} |\
        awk -v OFS='\\t' 'length($3)==1&&length($4)==1 {{ print $1,$2,$2+1,$5,"SNP";next }} {{ if ($3~/,/||$4~/,/) {{ print $1,$2,$2+1,$5,"MA" }} else {{ print $1,$2,$2+1,$5,"Indel" }} }}' |\
        bedtools intersect -wao -a - -b {input.annotation} |\
        awk -v OFS='\\t' -F '\\t' ' {{ ++A[$5][$9][$4] }} END {{ for (v in A) {{ for (c in A[v]) {{ for (t in A[v][c]) {{print v,c,t,A[v][c][t] }} }} }} }} ' > {output}
        '''

localrules: count_isec_overlaps
rule count_isec_overlaps:
    input:
        isec = expand(rules.bcftools_isec.output,chromosome=range(1,30),allow_missing=True),
    output:
        'vcfs/isec/{mode}.txt'
    shell:
        '''
        awk -v OFS='\\t' -F '\\t' ' {{ A[$1][$2][$3]+=$4 }} END {{ for (c in A) {{ for (k in A[c]) {{ for (t in A[c][k]) {{ print c,k,t,A[c][k][t] }} }} }} }}' {input} > {output}
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
