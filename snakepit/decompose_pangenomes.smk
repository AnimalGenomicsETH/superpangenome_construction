
rule vg_deconstruct:
    input:
        'graphs/{pangenome}/{chromosome}.gfa'
    output:
        'vcfs/{pangenome}/{chromosome}.raw.vcf'
    threads: 6
    resources:
        mem_mb= 5000,
        walltime= '4:00'
    params:
        ref_path = 'ref'
    shell:
        '''
        vg deconstruct -p {params.ref_path} -a -e -d 1 -t {threads} {input} |\
        sed 's/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t\\t/\\t.\\t/g;s/\\t$/\\t./g' {input}  |\
        awk '$0 !~ /ID=AT/{{print $0;next}}{{printf "%s\\n##INFO=<ID=CONFLICT,Number=1,Type=String,Description=Uncertain Paths>\\n", $0}}' \
        > {output}
        '''

rule vcfwave:
    input:
        rules.vg_deconstruct.output
    output:
        'vcfs/{pangenome}/{chromosome}.wavey.vcf'
    threads: 6
    resources:
        mem_mb= 5000,
        walltime= '24:00'
    params:
        skip_size = config.get('skip_size',100000)
    shell:
        '''
        vcfwave -t {threads} -L {params.skip_size} -k {input} > {output}
        '''

rule bcftools_norm:
    input:
        vcf = rules.vcfwave.output,
        ref = 'assemblies/{chromosome}/HER.fa'
    output:
        'vcfs/{pangenome}/{chromosome}.norm.vcf'
    shell:
        '''
        bcftools norm -m -any -f {input.ref} {input.vcf} |\
        bcftools norm -a |\
        bcftools norm -d none > {output}
        '''

rule bcftools_view:
    input:
        rules.bcftools_norm.output
    output:
        SV = 'vcfs/{pangenome}/{chromosome}.SV.vcf',
        small = 'vcfs/{pangenome}/{chromosome}.small.vcf'
    shell:
        '''
        bcftools view -i 'abs(ILEN)>=50' {input} > {output.SV}
        bcftools view -e 'abs(ILEN)>=50' {input} > {output.small}
        '''

asm_map = {'Angus':'asm5','Highland':'asm5','OBV':'asm5','BSW':'asm5','Pied':'asm5','Simmental':'asm5','Nellore':'asm10','Brahman':'asm10','Gaur':'asm20','Bison':'asm20','Yak':'asm20'}

rule minimap2_align:
    input:
        ref = 'assemblies/{chromosome}/HER.fa',
        query = 'assemblies/{chromosome}/{sample}.fa'
    output:
        temp('vcfs/assembly/{chromosome}.{sample}.paf')
    params:
        preset = lambda wildcards: config['pangenome_samples'][wildcards.sample]
    threads: 1
    resources:
        mem_mb= 10000,
        walltime= '4:00'
    shell:
        '''
        minimap2 -cx {params.preset} -t {threads} \
        --cs {input.ref} {input.query} |\
        sort -k6,6 -k8,8n > {output}
        '''

rule paftools_call:
    input:
        ref = 'assemblies/{chromosome}/HER.fa',
        paf = rules.minimap2_align.output
    output:
        temp('vcfs/assembly/{chromosome}.{sample}.raw.vcf')
    threads: 1
    resources:
        mem_mb = 2000,
        walltime = '1:00'
    shell:
        '''
        paftools.js call -f {input.ref} -s {wildcards.sample} {input.paf} > {output}
        '''

rule bcftools_merge:
    input:
        expand('vcfs/assembly/{{chromosome}}.{sample}.raw.vcf',sample=config['pangenome_samples'])
    output:
        temp('vcfs/assembly/{chromosome}.raw.vcf')
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools merge --threads {threads} -o {output} {input}
        '''
