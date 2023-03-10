
rule all:
    input:
        expand('working_area/{sample}_ASM/touch',sample=('BRA','ANG'))

rule bionano_denovo:
    input:
        bnx = 'working_area/{sample}.bnx',
        ref_cmap = 'working_area/reference1/ARS-UCD1.2_Btau5.0.1Y_DLE1_0kb_0labels.cmap'
    output:
        directory('working_area/{sample}_ASSEMBLY/contigs/exp_refineA'),
        'working_area/{sample}_ASSEMBLY/exp_optArguments.xml'
    params:
        _dir = 'working_area/{sample}_ASSEMBLY'
    envmodules:
        'gcc/4.8.5',
        'python/3.7.4',
        'perl/5.16.3',
        'r/4.1.3'
    threads: 16
    resources:
        mem_mb = 2500,
        walltime = '24h'
    shell:
        '''
        set +u
        . /cluster/apps/perl5/etc/bashrc

        mkdir -p {params._dir}/ClusterLogs

        ./runBNG.sh denovo -t /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/RefAligner/12432.12642rel -s /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25 -b {input.bnx} -T {threads} -j {threads} -l 50 -z 3000 -o {params} -r {input.ref_cmap} -a 16.2
        touch {output}
        '''

rule bionano_SV:
    input:
        query_cmaps = 'working_area/{sample}_ASSEMBLY/contigs/exp_refineA',
        args = 'working_area/{sample}_ASSEMBLY/exp_optArguments.xml',
        ref_cmap = 'working_area/reference1/ARS-UCD1.2_Btau5.0.1Y_DLE1_0kb_0labels.cmap'
    output:
        'working_area/{sample}_SV/touch'
    params:
        _dir = 'working_area/{sample}_SV'
    envmodules:
        'gcc/4.8.5',
        'python/3.7.4',
        'perl/5.16.3',
        'r/4.1.3'
    threads: 12
    resources:
        mem_mb = 5000,
        walltime = '24h'
    shell:
        '''
        python /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/Pipeline/20230127/runSV.py -t /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/RefAligner/12432.12642rel -r {input.ref_cmap} -q {input.query_cmaps} -o {params._dir} -p /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/Pipeline/20230127 -T {threads} -a {input.args}
        '''



rule bionano:
    input:
        bnx = 'working_area/{sample}.bnx',
        ref_cmap = 'working_area/reference1/ARS-UCD1.2_Btau5.0.1Y_DLE1_0kb_0labels.cmap',
        optArg = lambda wildcards: 'bionano/Angus_selected/Angus_nonhap_ES_sdb_m3_20180818_pipeline_results/exp_optArguments.xml' if wildcards.sample =='ANG' else 'bionano/Brahman_selected/Brahma_nonhap_ES_sdb_m3_20180818_pipeline_results/exp_optArguments.xml'
    output:
        'working_area/{sample}_ASM/touch'
        #directory('working_area/{sample}_ASSEMBLY/contigs/exp_refineA'),
        #'working_area/{sample}_ASSEMBLY/exp_optArguments.xml'
    params:
        _dir = 'working_area/{sample}_ASM'
    envmodules:
        'gcc/4.8.5',
        'python/3.7.4',
        'perl/5.16.3',
        'r/4.1.3'
    threads: 8
    resources:
        mem_mb = 2500,
        walltime = '24h'
    shell:
        '''
        set +u
        . /cluster/apps/perl5/etc/bashrc
        mkdir -p {params._dir}

        python /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/Pipeline/1.0/pipelineCL.py -R -d -U -T {threads} -j {threads} -J {threads} -TJ {threads} -Te {threads} -Tp {threads} -i 5 -t /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/RefAligner/12432.12642rel -l {params._dir} -a {input.optArg} -b {input.bnx} -r {input.ref_cmap} --species-reference non_human -finalmergeSV -guidedB -NoCheckFiles --autoRestart
        touch {output}
        '''

