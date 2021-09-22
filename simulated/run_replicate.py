
# config
norep = 1
nothresh = 1
listrep = list(range(1, norep + 1))
listthresh = list(range(1, nothresh + 1))
# listgr = ["pggb","cactus","minigraph"]
listgr = ["minigraph", "pggb"]
input_genome = "25_OBV.fa"
param_simfile = "paramsim_sv.tsv"

# directory
sifdir = "/cluster/work/pausch/danang/psd/bin/sif"
dirwork = "/cluster/work/pausch/danang/psd/scratch/sim/sim_rep"

rule all:
    input: "combine_result.tsv"

rule simulate_sv:
    input: input_genome
    output:
        "simulated_re_{rep}.fasta",
        "simulated_comb_{rep}.fasta",
        "simulated_{rep}.bed"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    params:
        simfile = param_simfile
    shell:
        """

        SURVIVOR simSV {input} {param_simfile} 0 0 simulated_{wildcards.rep}

        awk '$1 ~/>/ \
             {{$1=$1"_sim";print $0;next}} \
             {{print $0}}' simulated_{wildcards.rep}.fasta > {output[0]}

        cat {input} {output[0]} > {output[1]}

        """

rule construct_minigraph:
    input:
        backbone = input_genome,
        aug = rules.simulate_sv.output[0]
    output:
        graph = "minigraph_{rep}.gfa",
        vg = "minigraph_{rep}.vg",
        xg = "minigraph_{rep}.xg"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """
        minigraph -t {threads} -xggs {input.backbone} {input.aug} > {output.graph}

        vg convert -g {output.graph} -r 1 -v > {output.vg} 

        vg index -x {output.xg} {output.vg}

        """

rule deconstruct_minigraph:
    input: rules.construct_minigraph.output.xg
    output:
        "minigraph_decons_{rep}.vcf",
        "minigraph_decons_tidy_{rep}.vcf",
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "00:10"
    shell:
        """

        vg deconstruct -t {threads} -p 25_OBV {input} > {output[0]}

        awk  'BEGIN{{OFS="\\t"}} \
             $1 ~ /#/{{print $0; next}} \
             $4 >$5{{$8=$8";SVTYPE=DEL";print $0;next}} \
             $4 <$5{{$8=$8";SVTYPE=INS";print $0;next}}' {output[0]}  \
             > {output[1]}

        """

rule eval_minigraph:
    input:
        test = "minigraph_decons_tidy_{rep}.vcf",
        truth = "simulated_{rep}.bed"
    output:
        "result/comp_res_minigraph_rep{rep}_thresh{thresh}.tsv"
    threads: 5
    resources:
        mem_mb = 2000,
        walltime = "00:10"
    shell:
        """

        SURVIVOR eval {input.test} {input.truth} \
        {wildcards.thresh} mini_rep{wildcards.rep}_thresh{wildcards.thresh}

        """

rule construct_pggb:
    input:
        rules.simulate_sv.output[1]
    output:
        graph = touch("graph/pggb_{rep}/pggb_{rep}.finished")
    threads: 32
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    params:
        sifdir = sifdir,
        dirwork = dirwork
    shell:
        """

        singularity run --bind {params.dirwork} {params.sifdir}/pggb.sif \
        'pggb -i {input} -t {threads} -s 100000 -p 90 -n 10 \
        -v -L -S -m -o graph/pggb_{wildcards.rep}'

        """


rule convert_pggb:
    input: rules.construct_pggb.output
    output:
        vg = "graph/pggb_{rep}/pggb_{rep}.vg",
        xg = "graph/pggb_{rep}/pggb_{rep}.xg"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

        cd graph/pggb_{wildcards.rep}

        cp simulated_comb_{wildcards.rep}.fasta.*.smooth.gfa pggb_{wildcards.rep}.gfa

        vg convert -g pggb_{wildcards.rep}.gfa -x > pggb_{wildcards.rep}.xg
        vg convert -g pggb_{wildcards.rep}.gfa -v > pggb_{wildcards.rep}.vg

        """

rule deconstruct_graph:
    input:
        lambda wildcards: f"{wildcards.prog}_{wildcards.rep}.xg" if wildcards.prog == "minigraph" else f"graph/{wildcards.prog}_{wildcards.rep}/pggb_{wildcards.rep}.xg"
    output:
        "results/{prog}_decons_{rep}.vcf",
        "results/{prog}_decons_tidy_{rep}.vcf",
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "00:10"
    shell:
        """

        if [ {wildcards.prog} == "minigraph" ]; then 
             vg deconstruct -t {threads} -p 25_OBV {input} > {output[0]}
        else
             vg deconstruct -t {threads} -e -p 25_OBV {input} > {output[0]}
        fi

        awk  'BEGIN{{OFS="\\t"}} \
             $1 ~ /#/{{print $0; next}} \
             $4 >$5{{$8=$8";SVTYPE=DEL";print $0;next}} \
             $4 <$5{{$8=$8";SVTYPE=INS";print $0;next}}' {output[0]}  \
             > {output[1]}


        """

rule eval_graph:
    input:
        test = "results/{prog}_decons_tidy_{rep}.vcf",
        truth = "simulated_{rep}.bed"
    output:
        "results/comp_res_{prog}_rep{rep}_thresh{thresh}.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "00:10"
    shell:
        """

        SURVIVOR eval {input.test} {input.truth} \
        {wildcards.thresh} \
        {wildcards.prog}_rep{wildcards.rep}_thresh{wildcards.thresh}

        """

localrules: combine_results
rule combine_results:
    input:
        expand("results/comp_res_{prog}_rep{rep}_thresh{thresh}.tsv", prog=listgr, rep=listrep, thresh=listthresh)
    output:
        "combine_result.tsv"
    shell:
        """

        for file in {input}
        do
            prog = $(cut -f3 -d"_" <<< $file)
            rep = $(cut -f4 -d"_" <<< $file| sed 's/rep//')
            thresh = $(cut -f5 -d"_" <<< $file| sed 's/thresh//')
            awk -v prog=$prog -v rep=$rep -v thresh=$thresh \
            'END{{print prog,rep,thresh,$2,$3,$4,$5,$6,$7}}' >> {output}
        done
        """
