



rule vg_deconstruct:
    input:
        'graphs/{pangenome}/{chromosomes}.gfa'
    output:
        'vcfs/{pangenome}/{chromosome}.deconstruct.vcf'
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
        rules.vcfwave.output
    output:
        'vcfs/{pangenome}/{chromosome}.norm.vcf'
    shell:
        '''
        bcftools norm -m -any -f {params.reference} |\
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

