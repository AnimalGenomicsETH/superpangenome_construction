
rule all:
    input:
        expand('{sample}_SV_literal/{sample}.vcf',sample=('BRA','ANG'))
        #expand('working_area/{sample}_ASM3/touch',sample=('BRA','ANG'))

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
        'working_area/{sample}_ASM3/touch'
    params:
        _dir = 'working_area/{sample}_ASM3'
    envmodules:
        'gcc/4.8.5',
        'python/3.7.4',
        'perl/5.16.3',
        'r/4.1.3',
        'eth_proxy'
    threads: 24
    resources:
        mem_mb = 4500,
        walltime = '4h'
    shell:
        '''
        set +u
        . /cluster/apps/perl5/etc/bashrc
        mkdir -p {params._dir}

        python /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/Pipeline/1.0/pipelineCL.py -d -U -T {threads} -j {threads} -J {threads} -TJ {threads} -Te {threads} -Tp {threads} -Gsiz 2.7 -i 0 -A -N 1 -t /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/RefAligner/12432.12642rel -l {params._dir} -a {input.optArg} -b {input.bnx} -r {input.ref_cmap} --species-reference non_human -finalmergeSV -guidedB -NoCheckFiles --no-vaf --autoRestart
        touch {output}
        '''

rule direct:
    input:
        bnx = '/cluster/work/pausch/alex/Superpangenome/vcfs/bionano/working_area/{sample}.bnx'
    output:
        '{sample}_SV/touch'
    params:
        _dir = '/cluster/work/pausch/alex/Superpangenome/vcfs/bionano/{sample}_SV/{sample}'
    envmodules:
        'gcc/4.8.5',
        'python/3.7.4',
        'perl/5.16.3',
        'r/4.1.3',
        'eth_proxy'
    threads: 16
    resources:
        mem_mb = 3000,
        walltime = '24h'
    shell:
        '''
        set +u
        . /cluster/apps/perl5/etc/bashrc
        mkdir -p {params._dir}
        
        /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/RefAligner/12432.12642rel/RefAligner -ref /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/working_area/ANG_ASM3/ref/ARS-UCD1.2_Btau5.0.1Y_DLE1_0kb_0labels.cmap -maxthreads {threads} -f -output-veto-filter '(intervals.txt|\.map|\.maprate)$' -sv 5 -S 0.1 -T 1e-12 -indelMinConf 12 -A 12 -L 60 -hashgen 5 7 2.4 1.5 0.05 5.0 1 1 4 -hash -hashdelta 26 10 46 -hashoffset 1 -hashrange 0 -hashGC 300 -hashT2 1 -hashGrouped 5 7 -hashMultiMatch 30 -insertThreads {threads} -nosplit 2 -outlierBC -outlierExtend 12 120 -f -AlignRes 2.5 -ErrDist 2.0 -HSDrange 1.0 -Kmax 4 -MaxSE 0.5 -MinSF 0.05 -MinSR 0.01 -MinFP 0.05 -MinFN 0.005 -MinSD_Mult 0.5 -MultiMatches 2 1e-3 0.1 4 12.0 -MultiMatchesDelta 50.0 -MultiMatchesFilter 2 -endoutlier 3e-2 -outlier 3e-4 -CutFlip 3 0.2 1 10 1e-3 -9.9 1e-3 -9.9 3 -CutFlipMerge 1e-6 -InversionNormalPV 1e-10 -InversionOverlapFilter 0.7 15.0 -CutFlipSkip 20 -CutFlipFilter 0.50 0.25 0.99 -CutFlipBC 1 -MultiMatchesRev 1 -RefSplitStitch 1 -SmallDupFilter 0.2 0.8 0.2 1e-3 -CutFlipInvDup 5 0.5 -svDupSplit 3 30 3 50 10 -svInvDup_maxQryGapMult 0.8 -svInvDup_minRefOverlapFrac 0.5 -svDup_maxQryGapMult 0.8 -svDup_minRefOverlapSites 3 -MultiMatchesTotScore 2 12.0 -svInvDup_maxSpacing 120.0 -svDup_maxSpacing 200.0 -PVendoutlier -PVres 2 -RefSplit 1e-4 12 1e-5 -deltaX 12 -deltaY 12 -extend 1 -resEstimate -rres 0.9 -xmapUnique 5 -indelNoOverlap 3 -xmapchim 14 2000 -xmaplen -outlierType1 0 -outlierLambda 10 -svSideMinLen 20.0 -biaswt 0.2 -biaswtEnd 0.0 -biaswtOutlier 0.0 -svAlignConfByLen 0 -svMinConf 10 -finalsort-sitesdec -RAmem 3 1 -maxmem 248 -maxvirtmem 0 -FP 0.1 -FN 0.01 -o {params._dir} -i {input.bnx} -svSize 1
        touch {output}
        '''


rule direct2:
    input:
        bnx = '/cluster/work/pausch/alex/Superpangenome/vcfs/bionano/working_area/{sample}.bnx'
    output:
        smap = '{sample}_SV_literal/{sample}.smap',
        vcf = '{sample}_SV_literal/{sample}.raw.vcf'
    params:
        _dir = '/cluster/work/pausch/alex/Superpangenome/vcfs/bionano/{sample}_SV_literal/{sample}'
    envmodules:
        'gcc/4.8.5',
        'python/3.7.4',
        'perl/5.16.3',
        'r/4.1.3',
        'eth_proxy'
    threads: 24
    resources:
        mem_mb = 4000,
        walltime = '24h'
    shell:
        '''
        set +u
        . /cluster/apps/perl5/etc/bashrc
        mkdir -p {params._dir}

        /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/RefAligner/12432.12642rel/amd2/RefAligner -ref /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/working_area/ANG_ASM3/ref/ARS-UCD1.2_Btau5.0.1Y_DLE1_0kb_0labels.cmap -maxthreads {threads} -f -output-veto-filter '(intervals.txt|\.map|\.maprate)$' -sv 1 -FP 0.2 -FN 0.02 -sf 0.20 -sd 0.10 -mres 1e-3 -S 0.1 -T 1e-12 -indelMinConf 12 -A 8 -L 60 -hashgen 5 7 2.4 1.5 0.05 5.0 1 1 4 -hash -hashdelta 26 10 46 -hashoffset 1 -hashrange 0 -hashGC 300 -hashT2 1 -hashGrouped 5 7 -hashMultiMatch 30 -insertThreads 16 -nosplit 2 -outlierBC -outlierExtend 12 120 -f -AlignRes 2.5 -ErrDist 2.0 -HSDrange 1.0 -Kmax 4 -MaxSE 0.5 -MinSF 0.05 -MinSR 0.01 -MinFP 0.05 -MinFN 0.005 -MinSD_Mult 0.5 -MultiMatches 2 1e-3 0.1 4 12.0 -MultiMatchesDelta 50.0 -MultiMatchesFilter 2 -endoutlier 3e-2 -outlier 3e-4 -CutFlip 3 0.2 1 10 1e-3 -9.9 1e-3 -9.9 3 -CutFlipMerge 1e-6 -InversionNormalPV 1e-10 -InversionOverlapFilter 0.7 15.0 -CutFlipSkip 20 -CutFlipFilter 0.50 0.25 0.99 -CutFlipBC 1 -MultiMatchesRev 1 -RefSplitStitch 1 -SmallDupFilter 0.2 0.8 0.2 1e-3 -CutFlipInvDup 5 0.5 -svDupSplit 3 30 3 50 10 -svInvDup_maxQryGapMult 0.8 -svInvDup_minRefOverlapFrac 0.5 -svDup_maxQryGapMult 0.8 -svDup_minRefOverlapSites 3 -MultiMatchesTotScore 2 12.0 -svInvDup_maxSpacing 120.0 -svDup_maxSpacing 200.0 -PVendoutlier -PVres 2 -RefSplit 1e-4 12 1e-5 -deltaX 12 -deltaY 12 -extend 1 -resEstimate -rres 0.9 -xmapUnique 5 -indelNoOverlap 3 -xmapchim 14 2000 -xmaplen -outlierType1 0 -outlierLambda 10 -svSideMinLen 20.0 -biaswt 0.2 -biaswtEnd 0.0 -biaswtOutlier 0.0 -svAlignConfByLen 0 -svMinConf 1 -finalsort-sitesdec -RAmem 3 1 -maxmem 248 -maxvirtmem 0 -FP 0.33558 -FN 0.01399 -sf 0.14631 -sd -0.06769 -bpp 499.77 -res 1.971 -sr 0.0179 -o {params._dir} -i {input.bnx} -svSize 1

        python /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/tools/pipeline/Solve3.7_20221013_25/Pipeline/1.0/smap_to_vcf_v2.py -r /cluster/work/pausch/alex/Superpangenome/vcfs/bionano/working_area/ANG_ASM3/ref/ARS-UCD1.2_Btau5.0.1Y_DLE1_0kb_0labels.cmap -s {output.smap} -n {wildcards.sample} -o {output.vcf} -b False -i 1
        '''

rule clean_up:
    input:
        vcf = '{sample}_SV_literal/{sample}.raw.vcf'
    output:
        '{sample}_SV_literal/{sample}.vcf'
    shell:
        '''
        sed 's/chr//g' {input.vcf} > $TMPDIR/temp.vcf
        bcftools reheader -f /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai $TMPDIR/temp.vcf |\
        bcftools view -i '(SVTYPE=="INS"||SVTYPE=="DEL")&&ABS(INFO/SVLEN)<1000000' -c 1 -t 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 - > {output}
        '''
# awk -v OFS='\t' '!/#/&&$3<30&&$10>20 {print $3,int($7),int($8),$11,$10}' ANG_SV.indel | awk '($3-$2)>5000' | bedtools sort -i -
