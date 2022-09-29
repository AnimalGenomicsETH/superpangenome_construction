
pangenome_samples = config['pangenome_samples']

rule all:
    input:
        expand('graphs/{pangenome}/{chromosome}.gfa',pangenome=('minigraph','pggb','cactus'),chromosome=range(1,30))


rule minigraph_construct:
    input:
        order = 'sample_order.txt',
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa', sample=pangenome_samples)
    output: 
        temp('graphs/minigraph/{chromosome}.basic.gfa')
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '4:00'
    params:
        sample_order = lambda wildcards, input: ' '.join([l for l in open(input.order)]),
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
        assemblies = expand('assemblies/{{chromosome}}/{sample}_{{chromosome}}.fa', sample=pangenome_samples),
    output:
       gfa = "graph/pggb/{chromosome}.gfa",
    threads: 20
    resources:
        mem_mb = 10000,
        walltime = '24:00',
        disc_scratch = 30
    params:
        fasta = '$TMPDIR/all.fa.gz',
        divergence = config.get('divergence'),
        n_haplotypes = lambda wildcards, input: len(input.assemblies),
        sifdir = sifdir,
        dirwork = dirwork
    shell:
        '''
        cat {input} | bgzip -@ {threads} -c > {params.fasta}
        samtools faidx {params.fasta}
        
        singularity exec -B $TMPDIR {params.pggb} \
        'pggb -i {params.fasta} -t {threads} -s 100000 -p {params.divergence) -n {params.n_haplotypes} \
        -S -o graph'

        cp graph/*.smooth.gfa $LS_SUBCWD/{output.gfa}
        cp graph/*.smooth.og $LS_SUBCWD/{output.og}

        '''

localrules: cactus_seqfile
rule cactus_seqfile:
    input:
        assemblies = expand('assemblies/{{chromosome}}/{sample}_{{chromosome}}.fa.masked', sample=pangenome_samples),
        mash_distances = 'tree/distance.tri'
    output:
        temp('graph/cactus/{chromosome}.seqfile')
    threads: 1
    resources:
        mem_mb = 1000
    params:
        fasta_dir="assembly/{chromo}/{chromo}_rep"
    shell:
        '''
        ./construct_tree.R {input.distance_file} {params.fasta_dir} {output} fa.masked {wildcards.chromo}
        '''

rule cactus_construct:
    input:
        rules.cactus_seqfile.output
    output:
        temp('graph/cactus/{chromosome}.hal')
    threads: 20
    resources:
        mem_mb = 5000,
        walltime = "24:00"
    params:
        jobstore=jobstore + "/dump_{chromo}"
    shell:
        '''

        cactus --maxLocalJobs {threads} \
        {params.jobstore}  \
        {input} {output}

        '''


rule cactus_convert:
    input:
        rules.cactus_construct.output
    output:
        'graph/cactus/{chromosome}.gfa'
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
