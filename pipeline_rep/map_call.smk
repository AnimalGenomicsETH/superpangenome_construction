#!/usr/bin/env python
import os

norep = 6
listrep = list(range(1, norep + 1))
nothresh = 20
listthresh = list(range(1, nothresh + 1))
listlin = ["delly", "manta"]
listgraph = ["minigraph", "pggb"]
listall = listlin + listgraph
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
        bcf = "varcall/delly_{rep}/delly_{rep}_lin.bcf",
        vcf = "varcall/delly_{rep}/delly_{rep}_lin.vcf"
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
    output: "varcall/manta_{rep}/manta_{rep}_lin.vcf"
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

# localrules: eval_linear
# rule eval_linear:
    # input:
        # truth = "simulated/simulated_{rep}.bed",
        # test = "varcall/{prog}_{rep}/{prog}_{rep}_lin.vcf"
    # output: "results_map/res_{prog}_rep{rep}_thresh{thresh}.tsv"
    # shell:
        # """

        # SURVIVOR eval {input.test} {input.truth} \
        # {wildcards.thresh} \
        # results_map/{wildcards.prog}_rep{wildcards.rep}_thresh{wildcards.thresh} > {output}

        # """

rule modify_graph:
    input:
        graph = "graph/{prog}_{rep}/{prog}_{rep}.vg"
    output:
        graph = "graph/{prog}_{rep}/{prog}_{rep}_mod.vg",
        xg = "graph/{prog}_{rep}/{prog}_{rep}_mod.xg",
        gcsa = "graph/{prog}_{rep}/{prog}_{rep}_mod.gcsa"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """

        vg mod -X 256 {input.graph} > {output.graph}

        vg index -t {threads} -x {output.xg} -g {output.gcsa} {output.graph}

        """

rule map_graph:
    input:
        xg = "graph/{prog}_{rep}/{prog}_{rep}_mod.xg",
        gcsa = "graph/{prog}_{rep}/{prog}_{rep}_mod.gcsa",
        f1 = "sim_reads/simread_rep{rep}_1.fq",
        f2 = "sim_reads/simread_rep{rep}_2.fq"
    output:
        gam = "varcall/{prog}_{rep}/{prog}_{rep}.gam"
    threads: 32
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        vg map -t {threads} -x {input.xg} -g {input.gcsa} \
        -f {input.f1} {input.f2} > {output.gam}

        """

rule augment_variations:
    input:
        gam = "varcall/{prog}_{rep}/{prog}_{rep}.gam",
        graph = "graph/{prog}_{rep}/{prog}_{rep}_mod.vg"
    output:
        aug_vg = "varcall/{prog}_{rep}/{prog}_{rep}_aug.vg",
        aug_gam = "varcall/{prog}_{rep}/{prog}_{rep}_aug.gam",
        aug_xg = "varcall/{prog}_{rep}/{prog}_{rep}_aug.xg"
    threads: 10
    resources:
        mem_mb = 4000,
        walltime = "04:00"
    shell:
        """

        vg augment -t {threads} -m 5 -q 5 \
        {input.graph} {input.gam}  -A {output.aug_gam} > {output.aug_vg}

        vg index -t {threads} {output.aug_vg} -x {output.aug_xg}

        """


rule pack_gam:
    input:
        gam = "varcall/{prog}_{rep}/{prog}_{rep}_aug.gam",
        xg = "varcall/{prog}_{rep}/{prog}_{rep}_aug.xg"
    output:
        pack = "varcall/{prog}_{rep}/{prog}_{rep}_aug.pack"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """
        vg pack -t {threads} \
        -x {input.xg} -g {input.gam} -Q 5 -o {output.pack}
        """

rule find_snarls:
    input:
        xg = "varcall/{prog}_{rep}/{prog}_{rep}_aug.xg"
    output:
        snarls = "varcall/{prog}_{rep}/{prog}_{rep}_aug.snarl"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        vg snarls -t {threads} {input.xg} > {output.snarls}

        """

rule calling_graph:
    input:
        xg = "varcall/{prog}_{rep}/{prog}_{rep}_aug.xg",
        pack = "varcall/{prog}_{rep}/{prog}_{rep}_aug.pack",
        snarl = "varcall/{prog}_{rep}/{prog}_{rep}_aug.snarl",
    output:
        vcf = "varcall/{prog}_{rep}/{prog}_{rep}_graph.vcf"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        vg call {input.xg} -m 4,4 -t {threads} \
        -k {input.pack} -r {input.snarl} -p 25_OBV > {output.vcf}.temp

        awk '$1 ~ /#/ ||  $7 == "PASS"' {output.vcf}.temp |
        awk  'BEGIN{{OFS="\\t"}}$1 ~ /#/{{print $0; next}} \
        $4 >$5{{$8=$8";SVTYPE=DEL";print $0;next}} \
        $4 <$5{{$8=$8";SVTYPE=INS";print $0;next}}' > {output.vcf}


        """


def get_test_input(wildcards):
    if wildcards.prog in listlin:
        return "varcall/{prog}_{rep}/{prog}_{rep}_lin.vcf"
    else:
        return "varcall/{prog}_{rep}/{prog}_{rep}_graph.vcf"


localrules: eval_graph
rule eval_graph:
    input:
        truth = "simulated/simulated_{rep}.bed",
        test = get_test_input
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
        expand("results_map/res_{prog}_rep{rep}_thresh{thresh}.tsv", prog=listall, rep=listrep, thresh=listthresh)
    output:
        "results_map/combine_result_map.tsv"
    shell:
        """
        for file in {input}
        do
            prog=$(cut -f3 -d"_" <<< $file)
            rep=$(cut -f4 -d"_" <<< $file| sed 's/rep//')
            thresh=$(cut -f5 -d"_" <<< $file| sed -e 's/thresh//' -e 's/.tsv//')
            echo $prog $rep $thresh
            awk -v prog=$prog -v rep=$rep -v thresh=$thresh \
            'END{{print prog,rep,thresh,$2,$3,$4,$5,$6,$7}}' $file >> {output}
        done
        """
