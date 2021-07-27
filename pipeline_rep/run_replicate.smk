
# config
norep = 1
nothresh = 1
listrep = list(range(1, norep + 1))
listthresh = list(range(1, nothresh + 1))
listgr = ["pggb","cactus","minigraph"]
# listgr = ["minigraph", "pggb"]
#listgr = ["cactus"]
input_genome = "25_OBV.fa"
param_simfile = "paramsim_sv.tsv"

# directory
sifdir = "/cluster/work/pausch/danang/psd/bin/sif"
dirwork = "/cluster/work/pausch/danang/psd/scratch/sim/sim_rep2"
scrdir="/cluster/scratch/cdanang/scr_cac"

rule all:
    input: "combine_result.tsv"

rule simulate_sv:
    input: input_genome
    output:
        "simulated/simulated_re_{rep}.fasta",
        "simulated/simulated_comb_{rep}.fasta",
        "simulated/simulated_{rep}.bed"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    params:
        simfile = param_simfile
    shell:
        """

        SURVIVOR simSV {input} {params.simfile} 0 0 simulated/simulated_{wildcards.rep}

        awk '$1 ~/>/ \
             {{$1=$1"_sim";print $0;next}} \
             {{print $0}}' simulated/simulated_{wildcards.rep}.fasta > {output[0]}

        cat {input} {output[0]} > {output[1]}

        """

rule construct_minigraph:
    input:
        backbone = input_genome,
        aug = rules.simulate_sv.output[0]
    output:
        graph = "graph/minigraph_{rep}/minigraph_{rep}.gfa",
        vg = "graph/minigraph_{rep}/minigraph_{rep}.vg",
        xg = "graph/minigraph_{rep}/minigraph_{rep}.xg"
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
        -S -o graph/pggb_{wildcards.rep}'

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


rule cactus_graphmap:
        input: 
            graph="graph/minigraph_{rep}/minigraph_{rep}.gfa",
            fasta="simulated/simulated_re_{rep}.fasta"
        output: "graph/cactus_{rep}/cactus_{rep}.paf"
        threads: 10
        resources:
           mem_mb= 4000,
           walltime= "01:00"
        params:
            prefix="graph/cactus_{rep}",
            scrdir=scrdir +"/jobstore_map_{rep}"
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH

            echo -e "25_OBV\\t25_OBV.fa\\n25_sim\\t{input.fasta}" > {params.prefix}/seqfile_{wildcards.rep}.tsv 

            cactus-graphmap {params.scrdir}  \
            {params.prefix}/seqfile_{wildcards.rep}.tsv \
            {input.graph} {output} \
            --outputFasta {params.prefix}/cac_{wildcards.rep}.fa \
            --realTimeLogging

           """

rule cactus_align:
        input: "graph/cactus_{rep}/cactus_{rep}.paf"
        output: 
            gfa="graph/cactus_{rep}/cactus_{rep}.gfa.gz",
            hal="graph/cactus_{rep}/cactus_{rep}.hal",
            vg="graph/cactus_{rep}/cactus_{rep}.vg",
            xg="graph/cactus_{rep}/cactus_{rep}.xg"
        threads: 10
        resources:
           mem_mb= 4000,
           walltime= "01:00"
        params: 
            prefix="graph/cactus_{rep}",
            scrdir=scrdir +"/jobstore_align_{rep}"
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH

           cactus-align {params.scrdir} \
           {params.prefix}/seqfile_{wildcards.rep}.tsv \
           {input} {output.hal} \
           --pangenome --pafInput \
           --realTimeLogging \
           --outGFA --outVG --reference 25_OBV

           vg index -x {output.xg}  {output.vg} 

           """

rule deconstruct_graph:
    input: "graph/{prog}_{rep}/{prog}_{rep}.xg"
    output:
        "results/{prog}_decons_{rep}.vcf",
        "results/{prog}_decons_tidy_{rep}.vcf"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "00:10"
    shell:
        """

        if [[ "{wildcards.prog}" == "cactus" ]]; then 
            
             #vg deconstruct -t {threads} -e -p 25_OBV {input} > {output[0]}
             vg deconstruct -e -t {threads} \
             -p 25_OBV.25_OBV \
             -A "_MINIGRAPH" {input} > {output[0]} 

        else
             vg deconstruct -t {threads} -p 25_OBV {input} > {output[0]}
        fi

        #vg deconstruct -t {threads} -p 25_OBV {input} > {output[0]}


        if [[ "{wildcards.prog}" == "cactus" ]]; then

            awk  'BEGIN{{OFS="\\t"}} \
            $1 ~ /#/{{print $0; next}} \
            $4 > $5{{$1="25_OBV";$8=$8";SVTYPE=DEL";print $0;next}} \
            $4 < $5{{$1="25_OBV";$8=$8";SVTYPE=INS";print $0;next}}'  {output[0]}  \
            > {output[1]} 
        else

            awk  'BEGIN{{OFS="\\t"}} \
             $1 ~ /#/{{print $0; next}} \
             $4 > $5{{$8=$8";SVTYPE=DEL";print $0;next}} \
             $4 < $5{{$8=$8";SVTYPE=INS";print $0;next}}' {output[0]}  \
             > {output[1]}

        fi


        """

localrules: eval_graph
rule eval_graph:
    input:
        test = "results/{prog}_decons_tidy_{rep}.vcf",
        truth = "simulated/simulated_{rep}.bed"
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
        results/{wildcards.prog}_rep{wildcards.rep}_thresh{wildcards.thresh} > {output}

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
            prog=$(cut -f3 -d"_" <<< $file)
            rep=$(cut -f4 -d"_" <<< $file| sed 's/rep//')
            thresh=$(cut -f5 -d"_" <<< $file| sed 's/thresh//'| sed 's/.tsv//g')
            awk -v prog=$prog -v rep=$rep -v thresh=$thresh \
            'END{{print prog,rep,thresh,$2,$3,$4,$5,$6,$7}}' $file >> {output}
        done
        """