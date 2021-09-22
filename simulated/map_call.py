#!/usr/bin/env python
import os
from collections import defaultdict


norep = 1
listrep = list(range(1, norep + 1))
nothresh = 5
listthresh = list(range(1, nothresh + 1))
listgr = ["delly", "manta"]
rundir = os.getcwd()

rule all:
    input: "results_map/combine_result_map.tsv"

rule simulate_read:
    input:
        ref = "simulated/simulated_re_{rep}.fasta"
    output:
        f1 = "sim_reads/simread_rep{rep}_1.fq",
        f2 = "sim_reads/simread_rep{rep}_2.fq"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """
        # source $HOME/.bashrc
        # conda activate pan2

        mason_simulator -ir {input.ref} \
        -n 4500000  --num-threads {threads} -o {output.f1} -or {output.f2}

        """

rule map_linear:
    input:
        ref = "map_linear/25_OBV.fa",
        f1 = "sim_reads/simread_rep{rep}_1.fq",
        f2 = "sim_reads/simread_rep{rep}_2.fq"
    output:
        bam = "map_linear/lin_map_{rep}.bam",
        bai = "map_linear/lin_map_{rep}.bam.bai"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        bwa mem -t {threads} \
        {input.ref} {input.f1} {input.f2}| 
        samtools sort -@ {threads} -o {output.bam} -

        samtools index -@ {threads} {output.bam}

        """

rule call_delly:
    input:
        ref = "map_linear/25_OBV.fa",
        bam = "map_linear/lin_map_{rep}.bam"
    output:
        bcf = "varcall/delly_{rep}/delly_call_{rep}.bcf",
        vcf = "varcall/delly_{rep}/delly_call_{rep}.vcf"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        delly call -g {input.ref} {input.bam} -o {output.bcf} 
        bcftools view {output.bcf} > {output.vcf}

        """

rule call_manta:
    input:
        ref = "map_linear/25_OBV.fa",
        bam = "map_linear/lin_map_{rep}.bam"
    output: "varcall/manta_{rep}/manta_call_{rep}.vcf"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    params:
        rundir = rundir + "/varcall/manta_{rep}"
    shell:
        """

        # cd varcall/manta_{wildcards.rep}

        MANTA_INSTALL_PATH=/cluster/work/pausch/danang/psd/bin/manta-1.6.0.centos6_x86_64
        $MANTA_INSTALL_PATH/bin/configManta.py \
        --bam {input.bam} \
        --referenceFasta {input.ref} \
        --runDir {params.rundir}

        {params.rundir}/runWorkflow.py -j {threads}

        gunzip  {params.rundir}/results/variants/diploidSV.vcf.gz

        cp  {params.rundir}/results/variants/diploidSV.vcf {output}


        """

localrules: eval_linear
rule eval_linear:
    input:
        truth = "simulated/simulated_{rep}.bed",
        test = "varcall/{prog}_{rep}/{prog}_call_{rep}.vcf"
    output: "results_map/res_{prog}_rep{rep}_thresh{thresh}.tsv"
    shell:
        """

        SURVIVOR eval {input.test} {input.truth} \
        {wildcards.thresh} \
        results_map/{wildcards.prog}_rep{wildcards.rep}_thresh{wildcards.thresh} > {output}


        """
localrules: combine_results
rule combine_results:
    input:
        expand("results_map/res_{prog}_rep{rep}_thresh{thresh}.tsv", prog=listgr, rep=listrep, thresh=listthresh)
    output:
        "results_map/combine_result_map.tsv"
    shell:
        """
        for file in {input}
        do
            prog = $(cut -f2 -d"_" <<< $file)
            rep = $(cut -f3 -d"_" <<< $file| sed 's/rep//')
            thresh = $(cut -f4 -d"_" <<< $file| sed -e 's/thresh//' -e 's/.tsv//')
            awk -v prog=$prog -v rep=$rep -v thresh=$thresh \
            'END{{print prog,rep,thresh,$2,$3,$4,$5,$6,$7}}' >> {output}
        done
        """
