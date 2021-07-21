
norep = 6
listrep = list(range(2, norep + 1))
nothresh = 20
listthresh = list(range(1, nothresh + 1))
rundir = os.getcwd()
simtype = ["pacbio","ont"]


hmm_model_pacbio="/cluster/work/pausch/danang/psd/bin/pbsim2/data/P6C4.model"
hmm_model_ont="/cluster/work/pausch/danang/psd/bin/pbsim2/data/R95.model"
lin_ref = "25_OBV.fa"

rule all:
    input:
        "results_map/combine_result_call_lr_lin.tsv",


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
	   	  --length-min 1000 --length-max 100000 \
		  --hmm_model {params.hmm_model_ont}  --difference-ratio 23:31:46 \
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

localrules: eval_graph
rule eval_graph:
    input:
        truth = "simulated_sr/simulated_{rep}.bed",
        test = "call_linear/{sim}_{rep}_call.vcf"
    output: "results_map/res_lin_rep{rep}_thresh{thresh}_{sim}.tsv"
    shell:
        """

        SURVIVOR eval {input.test} {input.truth} \
        {wildcards.thresh} \
        results_map/lin_rep{wildcards.rep}_thresh{wildcards.thresh}_{wildcards.sim} > {output}


        """

localrules: combine_results
rule combine_results:
    input: expand("results_map/res_lin_rep{rep}_thresh{thresh}_{sim}.tsv",rep=listrep,thresh=listthresh,sim=simtype)
    output:
        "results_map/combine_result_call_lr_lin.tsv"
    shell:
        """
        for file in {input}
        do
            prog=$(cut -f3 -d"_" <<< $file)
            rep=$(cut -f4 -d"_" <<< $file| sed 's/rep//')
            thresh=$(cut -f5 -d"_" <<< $file| sed -e 's/thresh//')
			sim=$(cut -f6 -d"_" <<< $file|sed -e 's/.tsv//')
            echo $prog $rep $thresh
            awk -v prog=$prog -v rep=$rep -v thresh=$thresh -v sim=$sim \
            'END{{print prog,rep,thresh,sim,$2,$3,$4,$5,$6,$7}}' $file >> {output}
        done
        """
