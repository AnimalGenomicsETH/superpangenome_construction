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
        sed $'s/\\t\\t/\\t.\\t/g;s/\\t$/\\t./g' |\
        awk '$0 !~ /ID=AT/{{print $0;next}}{{printf "%s\\n##INFO=<ID=CONFLICT,Number=1,Type=String,Description=Uncertain Paths>\\n", $0}}' > {output}
        sed -i $'s/\t\t/\t.\t/g' {output}
        '''

rule vcfwave:
    input:
        rules.vg_deconstruct.output
    output:
        'vcfs/{pangenome}/{chromosome}.wavey.vcf'
    threads: 1
    resources:
        mem_mb= 15000,
        walltime= '120:00'
    params:
        skip_size = config.get('skip_size',0)
    shell:
        '''
        bcftools annotate -x INFO {input} |\
        vcfwave -t {threads} -L {params.skip_size} -k > {output}
        '''

rule bcftools_norm:
    input:
        vcf = rules.vcfwave.output,
        reference = expand('assemblies/{{chromosome}}/{reference_ID}.fa',reference_ID = get_reference_ID())
    output:
        'vcfs/{pangenome}/{chromosome}.norm.vcf'
    threads: 2
    resources:
        mem_mb = 5000,
        scratch = '5G'
    shell:
        '''
        bcftools norm --threads {threads} -m -any -f {input.reference} {input.vcf} |\
        bcftools norm --threads {threads} -d none |\
        bcftools sort -T $TMPDIR -o {output}
        #abcftools view -c 0 |\
        '''

#TODO probably want small vcf as gz for isec
rule bcftools_view:
    input:
        rules.bcftools_norm.output
    output:
        SV = 'vcfs/{pangenome}/{chromosome}.SV.vcf',
        small = multiext('vcfs/{pangenome}/{chromosome}.small.vcf.gz','','.tbi')
    threads: 1
    resources:
        mem_mb = 5000
    shell:
        '''
        bcftools view -i 'abs(ILEN)>=50' {input} > {output.SV}

        bcftools view -e 'abs(ILEN)>=50' {input} |\
        bcftools norm --threads {threads} -a |\
        bcftools norm --rm-dup exact |\
        bcftools view -c 1 |\
        bcftools norm -m +both |\
        bcftools sort -T $TMPDIR |\
        bcftools annotate -x INFO -o {output.small[0]}

        tabix -p vcf {output.small[0]}
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
        expand('vcfs/assembly/{{chromosome}}.{sample}.raw.vcf',sample=filter(lambda x: x != get_reference_ID(), pangenome_samples))
    output:
        temp('vcfs/assembly/{chromosome}.raw.vcf')
    threads: 2
    resources:
        mem_mb = 2500
    shell:
        '''
        bcftools merge --threads {threads} --no-index -o {output} {input}
        '''
