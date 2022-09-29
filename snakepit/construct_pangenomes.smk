


rule minigraph_construct:
    input:
        order = 'sample_order.txt',
        assemblies = expand('assemblies/{{chromosome}}/{sample}_{{chromosome}}.fa', sample=pangenome_samples)
    output: 'graphs/minigraph/{chromosome}.gfa'
    threads: 1
    resources:
        mem_mb = 10000,
        walltime = '4:00'
    params:
        sample_order = lambda wildcards: ' '.join([l for l in open(input.order)]),
        L = config.get('L',50),
        j = config.get('divergence',0.02)
    shell:
        '''
        minigraph -t {threads} -cxggs -j {params.j} -L {params.L} {params.sample_order} > {output}
        '''


rule pggb_construct:
    input:
        assemblies = expand('assemblies/{{chromosome}}/{sample}_{{chromosome}}.fa', sample=pangenome_samples)
    output:
       gfa = "graph/pggb/{chromosome}.gfa",
    threads: 20
    resources:
        mem_mb = 10000,
        walltime = "48:00",
        disc_scratch = 30
    params:
        fasta = '$TMPDIR/all.fa.gz',
        sifdir = sifdir,
        dirwork = dirwork
    shell:
        '''
        cat {input} | bgzip -@ {threads} -c > {params.fasta}
        samtools faidx {params.fasta}
        
        singularity exec -B $TMPDIR {params.pggb} \
        'pggb -i {params.fasta} -t {threads} -s 100000 -p 80 -n 10 \
        -S -o graph'

        cp graph/*.smooth.gfa $LS_SUBCWD/{output.gfa}
        cp graph/*.smooth.og $LS_SUBCWD/{output.og}

        '''


rule cactus_seqfile:
    input:
        assemblies = expand('assemblies/{{chromosome}}/{sample}_{{chromosome}}.fa.masked', sample=pangenome_samples),
        mash_distances = 'tree/distance.tri'
    output:
        temp('graph/cactus/{chromosome}.seqfile')
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "1:00"
    params:
        fasta_dir="assembly/{chromo}/{chromo}_rep"
    shell:
        '''
        ./construct_tree.R {input.distance_file} {params.fasta_dir} {output} fa.masked {wildcards.chromo}
        '''

rule cactus_construct:
    input:
        'graph/cactus/{chromosome}.seqfile'
    output:
        temp('graph/cactus/{chromosome}.hal')
    threads: 20
    resources:
        mem_mb = 5000,
        walltime = "24:00"
    params:
        jobstore=jobstore + "/dump_{chromo}"
    shell:
        """
        source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
        export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
        export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH

        cactus --maxLocalJobs {threads} \
        {params.jobstore}  \
        {input} {output}

        """


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
