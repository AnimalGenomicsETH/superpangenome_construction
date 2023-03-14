from pathlib import PurePath

rule all:
    input:
        expand('{sample}_SV/{sample}.vcf',sample=config['samples'])

rule bionano_fa2cmap:
    input:
        config['reference']
    output:
        f'reference/{PurePath(config["reference"]).with_suffix("").name}_DLE1_0kb_0labels.cmap'
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    envmodules:
        'gcc/4.8.5',
        'python/3.7.4',
        'perl/5.16.3',
        'r/4.1.3'
    resources:
        mem_mb = 10000
    shell:
        '''
        set +u
        . /cluster/apps/perl5/etc/bashrc
        perl {config[Solve]}HybridScaffold/20230127/scripts/fa2cmap_multi_color.pl -i {input} -o {params._dir} -e DLE1:CTTAAG 1
        '''

rule bionano_refaligner:
    input:
        bnx = '{sample}.bnx',
        cmap = rules.bionano_fa2cmap.output[0]
    output:
        smap = '{sample}_SV/{sample}.smap',
        vcf = '{sample}_SV/{sample}.raw.vcf'
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    envmodules:
        'gcc/4.8.5',
        'python/3.7.4',
        'perl/5.16.3',
        'r/4.1.3',
        'eth_proxy'
    threads: 24
    resources:
        mem_mb = 3000,
        walltime = '24h'
    shell:
        '''
        set +u
        . /cluster/apps/perl5/etc/bashrc
        mkdir -p {params._dir}

        {config[Solve]}RefAligner/12432.12642rel/RefAligner -ref {input.cmap} -maxthreads {threads} -f -output-veto-filter '(intervals.txt|\.map|\.maprate)$' -sv 1 -FP 0.2 -FN 0.02 -sf 0.20 -sd 0.10 -mres 1e-3 -S 0.1 -T 1e-12 -indelMinConf 12 -A 8 -L 60 -hashgen 5 7 2.4 1.5 0.05 5.0 1 1 4 -hash -hashdelta 26 10 46 -hashoffset 1 -hashrange 0 -hashGC 300 -hashT2 1 -hashGrouped 5 7 -hashMultiMatch 30 -insertThreads 16 -nosplit 2 -outlierBC -outlierExtend 12 120 -f -AlignRes 2.5 -ErrDist 2.0 -HSDrange 1.0 -Kmax 4 -MaxSE 0.5 -MinSF 0.05 -MinSR 0.01 -MinFP 0.05 -MinFN 0.005 -MinSD_Mult 0.5 -MultiMatches 2 1e-3 0.1 4 12.0 -MultiMatchesDelta 50.0 -MultiMatchesFilter 2 -endoutlier 3e-2 -outlier 3e-4 -CutFlip 3 0.2 1 10 1e-3 -9.9 1e-3 -9.9 3 -CutFlipMerge 1e-6 -InversionNormalPV 1e-10 -InversionOverlapFilter 0.7 15.0 -CutFlipSkip 20 -CutFlipFilter 0.50 0.25 0.99 -CutFlipBC 1 -MultiMatchesRev 1 -RefSplitStitch 1 -SmallDupFilter 0.2 0.8 0.2 1e-3 -CutFlipInvDup 5 0.5 -svDupSplit 3 30 3 50 10 -svInvDup_maxQryGapMult 0.8 -svInvDup_minRefOverlapFrac 0.5 -svDup_maxQryGapMult 0.8 -svDup_minRefOverlapSites 3 -MultiMatchesTotScore 2 12.0 -svInvDup_maxSpacing 120.0 -svDup_maxSpacing 200.0 -PVendoutlier -PVres 2 -RefSplit 1e-4 12 1e-5 -deltaX 12 -deltaY 12 -extend 1 -resEstimate -rres 0.9 -xmapUnique 5 -indelNoOverlap 3 -xmapchim 14 2000 -xmaplen -outlierType1 0 -outlierLambda 10 -svSideMinLen 20.0 -biaswt 0.2 -biaswtEnd 0.0 -biaswtOutlier 0.0 -svAlignConfByLen 0 -svMinConf 1 -finalsort-sitesdec -RAmem 3 1 -maxmem 248 -maxvirtmem 0 -FP 0.2 -FN 0.01 -sf 0.14631 -sd -0.06769 -bpp 499.77 -res 1.971 -sr 0.0179 -M 1 3 -o {params._dir} -i {input.bnx} -svSize 1

        python {config[Solve]}/Pipeline/20230127/smap_to_vcf_v2.py -r {input.cmap} -s {output.smap} -n {wildcards.sample} -o {params._dir}.raw -b False -i 1
        '''

rule clean_up_vcf:
    input:
        vcf = rules.bionano_refaligner.output['vcf']
    output:
        '{sample}_SV/{sample}.vcf'
    shell:
        '''
        sed 's/chr//g' {input.vcf} > $TMPDIR/temp.vcf
        bcftools reheader -f {config[reference]}.fai $TMPDIR/temp.vcf |\
        bcftools view -i '(SVTYPE=="INS"||SVTYPE=="DEL")&&ABS(INFO/SVLEN)<1000000' -c 1 -t 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 - > {output}
        '''
