
norep = 6
listrep = list(range(2, norep + 1))
nothresh = 20
listthresh = list(range(1, nothresh + 1))
rundir = os.getcwd()
simtype = ["pacbio","ont","hifi"]
listlin = ["lin"]
listgraph = ["minigraph", "pggb", "cactus"]
listall = listlin + listgraph

hmm_model_pacbio="/cluster/work/pausch/danang/psd/bin/pbsim2/data/P6C4.model"
hmm_model_ont="/cluster/work/pausch/danang/psd/bin/pbsim2/data/R95.model"
lin_ref = "25_OBV.fa"

rule all:
    input:
        "results_map/combine_result_call_lr_all.tsv",


rule simulate_long_read:
    input:
        genome="simulated_sr/simulated_re_{rep}.fasta"
    output:
        read="simulated_lr/{sim}_{rep}.fastq"
    threads: 10
    resources:
       mem_mb= 2000,
       walltime= "04:00"
    params:
       hmm_model_pacbio=hmm_model_pacbio,
       hmm_model_ont=hmm_model_ont
    shell:
       """
    
        if [[ {wildcards.sim} == "pacbio" ]]; then 

        echo "Simulate Pacbio data"
           pbsim2 --depth 20 \
             --length-min 1000 --length-max 25000 \
          --hmm_model {params.hmm_model_pacbio}  --prefix simulated_lr/{wildcards.sim}_{wildcards.rep} \
          {input.genome}

        elif [[ {wildcards.sim} == "ont" ]]; then
        echo "Simulate ONT data"

        pbsim2 --depth 20 \
          --length-min 5000 --length-max 100000 \
          --hmm_model {params.hmm_model_ont}  --difference-ratio 23:31:46 \
          --prefix simulated_lr/{wildcards.sim}_{wildcards.rep} \
          {input.genome}
        
        elif [[ {wildcards.sim} == "hifi" ]]; then
        echo "Simulate HiFi data"

        pbsim2 --depth 20 \
          --length-min 5000 --length-max 25000 \
          --hmm_model {params.hmm_model_pacbio} --accuracy-mean 0.98  \
          --prefix simulated_lr/{wildcards.sim}_{wildcards.rep} \
          {input.genome}

        fi


        mv simulated_lr/{wildcards.sim}_{wildcards.rep}_0001.fastq simulated_lr/{wildcards.sim}_{wildcards.rep}.fastq

       """

rule map_lr_linear:
    input:
        read="simulated_lr/{sim}_{rep}.fastq"
    output: "call_linear/{sim}_{rep}_map.bam"
    threads: 32
    resources:
       mem_mb= 2000,
       walltime= "04:00"
    params:
        lin_ref=lin_ref
    shell:
       """

       ngmlr -t {threads} -r {params.lin_ref} -q {input.read} -o temp_{wildcards.sim}_{wildcards.rep}.bam

       samtools sort -@ {threads} temp_{wildcards.sim}_{wildcards.rep}.bam | 
       samtools view -hb > {output}

       rm temp_{wildcards.sim}_{wildcards.rep}.bam

       """

rule call_sniffles:
    input:
        "call_linear/{sim}_{rep}_map.bam"
    output:
        "call_linear/{sim}_{rep}_call.vcf"
    threads: 32
    resources:
       mem_mb= 2000,
       walltime= "01:00"
    shell:
       """

       sniffles -t {threads} -m {input} -v {output}

       """

rule make_compatible_gfa:
    input: 
        lambda wildcards: "graph/{prog}_{rep}/{prog}_{rep}.gfa.gz" if wildcards.prog == "cactus" else "graph/{prog}_{rep}/{prog}_{rep}.gfa"
    output: 
         gfa = "graph/{prog}_{rep}/{prog}_comp_{rep}.gfa",
         xg = "graph/{prog}_{rep}/{prog}_comp_{rep}.xg"
    threads: 10
    resources:
       mem_mb= 2000,
       walltime= "01:00"
    shell:
       """
       
       if [[ {wildcards.prog} == "minigraph" ]]; then

            vg convert -g {input} -r 0 -f > {output.gfa}
            vg convert -g {output.gfa} -x > {output.xg}

       elif [[ {wildcards.prog} == "cactus" ]]; then

            zless {input} |
            awk '$1 !~/P/{{print;next}} \
            $0 ~ /OBV/ {{split($2,arr,"\.");print $1, arr[1],$3,$4}}' OFS="\\t" > {output.gfa}
            vg convert -g {output.gfa} -x > {output.xg}
        
       elif [[ {wildcards.prog} == "pggb" ]]; then

             ln --relative -s {input} {output.gfa}
             vg convert -g {output.gfa} -x > {output.xg}
       fi

       """

rule map_graph:
    input: 
        graph = "graph/{prog}_{rep}/{prog}_comp_{rep}.gfa",
        read = "simulated_lr/{sim}_{rep}.fastq"
    output: 
        gam = "call_graph/{prog}_{rep}/{prog}_{rep}_{sim}.gam"
    threads: 32
    resources:
       mem_mb= 2000,
       walltime= "04:00"
    shell:
       """

       module load gcc/6.3.0

       GraphAligner -t {threads} -x vg \
       --graph {input.graph} -f {input.read} -a {output.gam}
       
       """

rule pack_graph:
    input:
        xg = "graph/{prog}_{rep}/{prog}_comp_{rep}.xg",
        gam = "call_graph/{prog}_{rep}/{prog}_{rep}_{sim}.gam"
    output: "call_graph/{prog}_{rep}/{prog}_{rep}_{sim}.pack"
    threads: 32
    resources:
       mem_mb= 2000,
       walltime= "01:00"
    shell:
       """

           vg pack -t {threads} -x {input.xg} \
           -g {input.gam} \
           -o {output}

       """

rule call_graph:
    input:
        xg = "graph/{prog}_{rep}/{prog}_comp_{rep}.xg",
        pack = "call_graph/{prog}_{rep}/{prog}_{rep}_{sim}.pack"
    output: "call_graph/{prog}_{rep}/{prog}_{rep}_{sim}.vcf"
    threads: 32
    resources:
       mem_mb= 2000,
       walltime= "01:00"
    shell:
       """
        vg call -m 4,4 -t {threads} -k {input.pack} {input.xg} > {output}.temp

        awk '$1 ~ /#/ ||  $7 == "PASS" {{print}}' {output}.temp|
        awk  'BEGIN{{OFS="\t"}}$1 ~ /#/{{print $0; next}} 
            $4 >$5{{$8=$8";SVTYPE=DEL";print $0;next}}
            $4<$5{{$8=$8";SVTYPE=INS";print $0;next}}' > {output}

        rm {output}.temp


       """

def get_eval_test_data(wildcards):
    if wildcards.prog in listgraph:
        return "call_graph/{prog}_{rep}/{prog}_{rep}_{sim}.vcf"
    else:
        return "call_linear/{sim}_{rep}_call.vcf"


localrules: eval_graph
rule eval_graph:
    input:
        truth = "simulated_sr/simulated_{rep}.bed",
        test = get_eval_test_data
    output: "results_map/res_{prog}_rep{rep}_thresh{thresh}_{sim}.tsv"
    shell:
        """

        SURVIVOR eval {input.test} {input.truth} \
        {wildcards.thresh} \
        results_map/{wildcards.prog}_rep{wildcards.rep}_thresh{wildcards.thresh}_{wildcards.sim} > {output}


        """

localrules: combine_results
rule combine_results:
    input: expand("results_map/res_{prog}_rep{rep}_thresh{thresh}_{sim}.tsv",prog=listall, rep=listrep,thresh=listthresh,sim=simtype)
    output:
        "results_map/combine_result_call_lr_all.tsv"
    shell:
        """
        for file in {input}
        do
            prog=$(cut -f3 -d"_" <<< $file)
            rep=$(cut -f4 -d"_" <<< $file| sed 's/rep//')
            thresh=$(cut -f5 -d"_" <<< $file| sed -e 's/thresh//')
            sim=$(cut -f6 -d"_" <<< $file|sed -e 's/.tsv//')
            echo $prog $rep $thresh $sim
            awk -v prog=$prog -v rep=$rep -v thresh=$thresh -v sim=$sim \
            'END{{print prog,rep,thresh,sim,$2,$3,$4,$5,$6,$7}}' $file >> {output}
        done
        """
