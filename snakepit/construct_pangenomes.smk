
pangenome_samples = config['pangenome_samples']

wildcard_constraints:
    chromosome = r'\d+',
    pangenome = r'pggb|cactus|minigraph|assembly'

include: 'decompose_pangenomes.smk'

rule all:
    input:
        expand('graphs/{pangenome}/{chromosome}.gfa',pangenome=('minigraph','pggb','cactus'),chromosome=range(1,30)),
        expand('vcfs/{pangenome}/{chromosome}.SV.vcf',pangenome=('minigraph','pggb','cactus','assembly'),chromosome=range(1,30))

##normalise direction of contigs

#rule download_fasta:
#    output:
#        expand('raw_assemblies/{sample}.fasta',sample=config['pangenome_samples'])
#    shell:
#        '''
#        ./prep_assemblies.sh
#        '''

rule split_fasta:
    input:
        'raw_assemblies/{sample}.fasta'
    output:
        'assemblies/{chromosome}/{sample}.fa'
    shell:
        '''
        echo "{wildcards.chromosome}" | seqtk subseq {input} - |\
        sed 's/>.*/>{wildcards.sample}/' > {output}
        '''

rule repeatmasker_soft:
    input:
        'assemblies/{chromosome}/{sample}.fa'
    output:
        temp('assemblies/{chromosome}/{sample}.fa.masked')
    threads: 8
    resources:
        mem_mb = 1500,
        walltime = '4:00'
    shell:
        '''
        RepeatMasker -pa $(({threads}/2)) -no_is -qq -xsmall \
        -lib {config[repeat_library]} {input}
        '''

rule mash_triangle:
    input:
        lambda wildcards: expand('raw_assemblies/{sample}.fasta',sample=config[wildcards.sample_set])
    output:
        'tree/{sample_set}.lower_triangle.txt'
    threads: 6
    resources:
        mem_mb = 2500,
        walltime = '4:00'
    shell:
        '''
        mash triangle -s 10000 -k 25 -p {threads} {input} | awk 'NR>1' {output}
        '''        

rule minigraph_construct:
    input:
        mash_distances = expand(rules.mash_triangle.output,sample_set='pangenome_samples'),
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa', sample=pangenome_samples)
    output: 
        temp('graphs/minigraph/{chromosome}.basic.gfa')
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '4:00'
    params:
        sample_order = lambda wildcards, input: ' '.join([l for l in open(input.mash_distances[0])]),
        L = config.get('L',50),
        j = config.get('divergence',0.02)
    shell:
        '''
        minigraph -t {threads} -cxggs -j {params.j} -L {params.L} {params.sample_order} > {output}
        '''

rule minigraph_call:
    input:
        gfa = rules.minigraph_construct.output,
        sample = 'assemblies/{chromosome}/{sample}.fa'
    output:
        'graphs/minigraph/{chromosome}.{sample}.bed'
    params:
        L = config.get('L',50),
        j = config.get('divergence',0.02)
    shell:
        '''
        minigraph -t {threads} -cxasm --call -j {params.j} -L {params.L} {input.gfa} {input.sample} > {output}
        '''

rule minigraph_path:
    input:
        paths = expand('graphs/minigraph/{{chromosome}}.{sample}.bed',sample=pangenome_samples),
        gfa = 'graphs/minigraph/{chromosome}.basic.gfa'
    output:
        'graphs/minigraph/{chromosome}.gfa'
    shell:
        '''
        add p lines
        '''

rule pggb_construct:
    input:
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa', sample=pangenome_samples),
    output:
       gfa = 'graphs/pggb/{chromosome}.gfa',
    threads: 20
    resources:
        mem_mb = 10000,
        walltime = '24:00',
        disc_scratch = 30
    params:
        pggb = '',
        fasta = '$TMPDIR/all.fa.gz',
        divergence = config.get('divergence'),
        n_haplotypes = lambda wildcards, input: len(input.assemblies),
    shell:
        '''
        cat {input} | bgzip -@ {threads} -c > {params.fasta}
        samtools faidx {params.fasta}
        
        singularity exec -B $TMPDIR {params.pggb} \
        pggb -i {params.fasta} -t {threads} -s 100000 -p {params.divergence} -n {params.n_haplotypes} \
        -S -o graph

        cp graph/*.smooth.gfa {output.gfa}

        '''

localrules: cactus_seqfile
rule cactus_seqfile:
    input:
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa.masked', sample=pangenome_samples),
        mash_distances = expand(rules.mash_triangle.output,sample_set='pangenome_samples')
    output:
        temp('graphs/cactus/{chromosome}.seqfile')
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        '''
        make_tree
        '''

rule cactus_construct:
    input:
        rules.cactus_seqfile.output
    output:
        temp('graphs/cactus/{chromosome}.hal')
    threads: 20
    resources:
        mem_mb = 5000,
        walltime = '24:00'
    shell:
        '''

        cactus --maxLocalJobs {threads} \
        {input} {output}

        '''

rule cactus_convert:
    input:
        rules.cactus_construct.output
    output:
        'graphs/cactus/{chromosome}.gfa'
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '4:00'
    shell:
        '''
        hal2vg --inMemory {input} |\
        vg convert -t {threads} -f |\
        awk '$1 !~ /P/{{print;next}} $2 !~ /Anc/{{print}}' > {output}
        '''
