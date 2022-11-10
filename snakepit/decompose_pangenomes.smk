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
        ref_path = get_reference_ID()
    shell:
        '''
        vg deconstruct -p {params.ref_path} -a -e -d 1 -t {threads} {input} |\
        awk '$0 !~ /ID=AT/{{print $0;next}}{{printf "%s\\n##INFO=<ID=CONFLICT,Number=1,Type=String,Description=Uncertain Paths>\\n", $0}}' > {output}
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
        reference = expand('assemblies/{{chromosome}}/{reference_ID}.fa',reference_ID = get_reference_ID())
    output:
        'vcfs/{pangenome}/{chromosome}.norm.vcf'
    shell:
        '''
        bcftools norm -m -any -f {input.reference} {input.vcf} |\
        bcftools norm -a |\
        bcftools norm -d none > {output}
        '''

#TODO probably want small vcf as gz for isec
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

rule minimap2_align:
    input:
        reference = expand('assemblies/{{chromosome}}/{reference_ID}.fa',reference_ID = get_reference_ID()),
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
        --cs {input.reference} {input.query} |\
        sort -k6,6 -k8,8n > {output}
        '''

rule paftools_call:
    input:
        reference = expand('assemblies/{{chromosome}}/{reference_ID}.fa',reference_ID = get_reference_ID()),
        paf = rules.minimap2_align.output
    output:
        temp('vcfs/assembly/{chromosome}.{sample}.raw.vcf')
    threads: 1
    resources:
        mem_mb = 2000,
        walltime = '1:00'
    shell:
        '''
        paftools.js call -f {input.reference} -s {wildcards.sample} {input.paf} > {output}
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
