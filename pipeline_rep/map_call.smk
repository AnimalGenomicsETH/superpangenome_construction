#!/usr/bin/env python
import os

norep = 6
listrep = list(range(2, norep + 1)) 
#add non-var
listrep.append(23)
nothresh = 20
listthresh = list(range(1, nothresh + 1))
listlin = ["delly", "manta"]
listgraph = ["minigraph", "pggb", "cactus","vglin"]
listall = listlin + listgraph
rundir = os.getcwd()

rule all:
    input:
        "results_map/combine_result_mapping.tsv",
        "results_map/combine_result_call.tsv"

rule simulate_read:
    input:
        ref = "simulated/simulated_re_{rep}.fasta"
    output:
        f1 = "sim_reads/simread_rep{rep}_1.fq",
        f2 = "sim_reads/simread_rep{rep}_2.fq",
        #tpos = "sim_reads/simread_rep{rep}_pos.sam"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """
        # source $HOME/.bashrc
        # conda activate pan2

        mason_simulator -ir {input.ref} \
        -n 4500000  --num-threads {threads} \
        --illumina-read-length 150 \
        -o {output.f1} -or {output.f2} -oa {output.tpos}

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
    params:
        prefix = "graph/{prog}_{rep}"
    shell:
        """

        if [[ {wildcards.prog} == "cactus" ]]; then

            zless {params.prefix}/cactus_{wildcards.rep}.gfa.gz |
            awk '$1 !~/P/{{print;next}} \
            $0 ~ /OBV/ {{split($2,arr,"\.");print $1, arr[1],$3,$4}}' \
            OFS="\\t" > {params.prefix}/temp.gfa 

            vg convert -g {params.prefix}/temp.gfa -v > {params.prefix}/temp.vg

            vg mod -X 256 {params.prefix}/temp.vg > {output.graph}

            vg index -t {threads} -x {output.xg} -g {output.gcsa} {output.graph}
        else  
            vg mod -X 256 {input.graph} > {output.graph}

            vg index -t {threads} -x {output.xg} -g {output.gcsa} {output.graph}
        fi

        """

#construct linear graph
rule construct_graph_linear:
        input: ref = "map_linear/25_OBV.fa"
        output: "graph/vglin_{rep}/vglin_{rep}_mod.vg"
        threads: 10
        resources:
           mem_mb= 2000,
           walltime= "00:30"
        shell:
           """

           vg construct -t {threads} -m 256 -r {input.ref} > {output}

           """

rule index_graph_linear:
        input: 
            graph = "graph/vglin_{rep}/vglin_{rep}_mod.vg"
        output: 
            xg = "graph/vglin_{rep}/vglin_{rep}_mod.xg",
            gcsa = "graph/vglin_{rep}/vglin_{rep}_mod.gcsa"
        threads: 10
        resources:
           mem_mb= 5000,
           walltime= "04:00"
        shell:
           """

           vg index -t {threads} -x {output.xg} -g {output.gcsa} {input.graph}

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
        -f {input.f1} -f {input.f2} > {output.gam}

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
        "results_map/combine_result_call.tsv"
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

# mapping statistics
rule stat_map_linear:
    input:
        bam = "map_linear/lin_map_{rep}.bam"
    output: "map_linear/map_stat_{rep}.tsv"
    threads: 1
    resources:
        mem_mb = 1000,
        walltime = "00:10"
    shell:
        """

        # samtools stats --threads {threads} {input.bam} > {output}

        ./bam_stat.py -i {input.bam} -o {output}

        """

rule stat_map_graph:
    input:
        gam = "varcall/{prog}_{rep}/{prog}_{rep}.gam",
        xg = "graph/{prog}_{rep}/{prog}_{rep}_mod.xg"
    output: "varcall/{prog}_{rep}/{prog}_{rep}_stat.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    params:
        prefix = "varcall/{prog}_{rep}"
    shell:
        """

        vg stats -a {input.gam} > {params.prefix}/temp1.tsv
        tmp1={params.prefix}/temp1.tsv

        vg convert -t 10 --gam-to-gaf {input.gam} \
        {input.xg} > {params.prefix}/tmp.gaf

        totalign=$(grep "Total alignments" $tmp1|cut -f3 -d" ")
        prmap=$(grep "Total aligned" $tmp1|cut -f3 -d" ")
        unmap=$(( totalign-prmap ))
        mq0=$(awk '$12==0' {params.prefix}/tmp.gaf  |wc -l)
        mq10=$(awk '$12>=10' {params.prefix}/tmp.gaf  |wc -l)
        mq60=$(awk '$12==60' {params.prefix}/tmp.gaf  |wc -l)
        al99=$(grep "perfect" $tmp1|cut -f3 -d" ")
        alp=$(grep "perfect" $tmp1|cut -f3 -d" ")
        clip=$(grep "Softclips" $tmp1|cut -f5 -d" ")
        pp=$(grep "properly" $tmp1|cut -f4 -d" ")


        echo graph {wildcards.prog} {wildcards.rep} \
             $totalign $unmap $mq0 $mq10 $mq60 \
             $al99 $alp $clip $pp > {output}


        rm {params.prefix}/tmp.gaf $tmp1

        """


def get_map_infile(listgraph, listrep):
    # add linear
    all_input = [f"map_linear/map_stat_{x}.tsv" for x in listrep]
    # add graph
    for comp in listgraph:
        all_input.extend(f"varcall/{comp}_{x}/{comp}_{x}_stat.tsv" for x in listrep)
    return all_input


localrules: combine_map
rule combine_map:
    input: get_map_infile(listgraph, listrep)
    output: "results_map/combine_result_mapping.tsv"
    shell:
        """

        cat {input} > {output}

        """
