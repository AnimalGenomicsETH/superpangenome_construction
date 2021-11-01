#!/usr/bin/env python

dirwork="/cluster/work/pausch/danang/psd/scratch/real_yak2"
#chromo_list=[2,3] + list(range(5,30))
chromo_list=[25]
assemb_list = ["UCD","Angus","Highland","OBV","Brahman","Yak","BSW","Pied","Gaur","Nellore"]
#assemb_list = ["UCD","Angus"]
ref_genome = "UCD"
rep_library="BosTau9_repeat_library.fasta"
#temporary files for cactus
jobstore="/cluster/scratch/cdanang/cactus_dump"
#directory for pggb
dirwork = "/cluster/work/pausch/danang/psd/scratch/real_all"
sifdir = "/cluster/work/pausch/danang/psd/bin/sif"


def sample_def(sample_file):
    allcomb=dict()
    with open(sample_file) as infile:
        for line in infile:
            token=line.strip().split()
            allcomb[token[0]]=token[1].split(",")
    return allcomb

allcomb=sample_def("sample_design.tsv")
grname=list(allcomb.keys())
allassemb=[]

for value in allcomb.values():
    for assemb in value:
        if assemb not in allassemb:
            allassemb.append(assemb)

print(allassemb)
print(grname)


rule all:
    input:
        expand("assembly/{chromo}/{chromo}_rep/{breed}_{chromo}.fa.masked",breed=allassemb,chromo=chromo_list),
        #expand("tree/combined_sketch_{grtype}.msh",grtype=grname),
        #expand("graph/cactus/cactus_{chromo}_{grtype}_seqfile.tsv",grtype=grname,chromo=chromo_list)
        #expand("graph/cactus/cactus_{chromo}_{grtype}_seqfile.tsv",grtype=grname,chromo=chromo_list)
        ##expand("graph/cactus/cactus_{chromo}_{grtype}.hal",grtype=grname,chromo=chromo_list),
	#expand("graph/minigraph/minigraph_{chromo}_{grtype}.gfa",grtype=grname,chromo=chromo_list),
	expand("graph/pggb/pggb_{chromo}_{grtype}/pggb_{chromo}_{grtype}.gfa",grtype=grname,chromo=chromo_list)


    
rule cactus_masking_chromosome:
    input:"assembly/{chromo}/{breed}_{chromo}.fa"
    output:"assembly/{chromo}/{chromo}_rep/{breed}_{chromo}.fa.masked"
    threads: 10
    resources:
        mem_mb= 1000,
        walltime= "08:00"
    params:
        rep_library=rep_library
    shell:
        """

        RepeatMasker -dir assembly/{wildcards.chromo}/{wildcards.chromo}_rep -no_is \
                -qq -xsmall -lib {params.rep_library} {input}

        """

rule sketch_assembly:
    input: "assembly/{breed}_aut.fa"
    output: "tree/{breed}.fa.msh"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        mash sketch -p {threads} -o {output} {input}

        """



rule combined_sketch:
    #input: expand("tree/{breed}.fa.msh", breed=assemb_list)
    input: lambda wildcards: [f"tree/{breed}.fa.msh" for breed in allcomb[wildcards.grtype]]
    output: "tree/combined_sketch_{grtype}.msh"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        mash paste {output} {input}

        """

rule estimate_distance:
    input: "tree/combined_sketch_{grtype}.msh"
    output: "tree/combined_distance_{grtype}.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

        mash dist {input} {input} > {output}

        """


rule cactus_create_seqfile:
    input:
        assemb=lambda wildcards: [f"tree/{breed}.fa.msh" for breed in allcomb[wildcards.grtype]],
        distance_file=rules.estimate_distance.output
    output:"graph/cactus/cactus_{chromo}_{grtype}_seqfile.tsv"
    threads: 10
    resources:
        mem_mb=1000,
        walltime="01:00"
    params:
        fasta_dir="assembly/{chromo}/{chromo}_rep"
    shell:
        """

        ./construct_tree.R {input.distance_file} {params.fasta_dir} {output} fa.masked {wildcards.chromo}

        """


rule cactus_align:
    input:"graph/cactus/cactus_{chromo}_{grtype}_seqfile.tsv"
    output:"graph/cactus/cactus_{chromo}_{grtype}.hal"
    threads: 20
    resources:
        mem_mb=5000,
        walltime="24:00"
    params:
        jobstore=jobstore + "/dump_{chromo}_{grtype}"
    shell:
        """


        source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
        export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
        export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH

        cactus --maxLocalJobs {threads} \
        {params.jobstore}  \
        {input} {output}

        """
        


### minigraph

rule determine_order:
    input: "tree/combined_distance_{grtype}.tsv"
    output: "tree/order_breed/breed_order_{grtype}.tsv"
    threads: 10
    resources:
        mem_mb=1000,
        walltime="01:00"
    params:
        ref = ref_genome
    run:
        import pandas as pd

        allin = pd.read_csv(input[0], sep="\t", header=None)
        allin.columns = ["anim1", "anim2", "dist", "comp1", "comp2"]
        allin["anim1"] = allin.anim1.str.extract(r"([A-Z][^_]+)")
        allin["anim2"] = allin.anim2.str.extract(r"([A-Z][^_]+)")

        allsel = allin.loc[allin.anim1 == params.ref]

        allsort = allsel.sort_values(by="dist")
        outfile = open(output[0], "w")
        print(*list(allsort.anim2), file=outfile)


def get_order_list(grtype, chromo):
    print(grtype,chromo)
    order_file=f"tree/order_breed/breed_order_{grtype}.tsv"
    if not os.path.exists(order_file): return ""
    infile = open(order_file)
    print(f"assembly/{chromo}/{br}_{chromo}.fa" for br in infile.readlines()[0].strip().split(" "))
    return [f"assembly/{chromo}/{br}_{chromo}.fa" for br in infile.readlines()[0].strip().split(" ")]
    infile.close()


rule run_minigraph:
    input:
        order_file = "tree/order_breed/breed_order_{grtype}.tsv",
        assembly = lambda wildcards: [f"assembly/{wildcards.chromo}/{breed}_{wildcards.chromo}.fa" for breed in allcomb[wildcards.grtype]]
    output: "graph/minigraph/minigraph_{chromo}_{grtype}.gfa"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    params: order_list = get_order_list("{grtype}", "{chromo}")
    shell:
        """
	echo {wildcards.grtype}	
        echo {params.order_list}

        minigraph -t {threads} -xggs {params.order_list} > {output}

        """


#pggb 

rule combine_genome:
    input: lambda wildcards: [f"assembly/{wildcards.chromo}/{breed}_{wildcards.chromo}.fa" for breed in allcomb[wildcards.grtype]]
    output: "comb_chromo/{chromo}_{grtype}_comb.fa"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "00:10"
    shell:
        """

        cat {input} > {output}

        """

rule construct_pggb:
    input:
        "comb_chromo/{chromo}_{grtype}_comb.fa"
    output:
        "graph/pggb/pggb_{chromo}_{grtype}/pggb_{chromo}_{grtype}.gfa"
    threads: 32
    resources:
        mem_mb = 2000,
        walltime = "24:00"
    params:
        sifdir = sifdir,
        dirwork = dirwork
    shell:
        """

        singularity run --bind {params.dirwork} {params.sifdir}/pggb.sif \
        'pggb -i {input} -t {threads} -s 100000 -p 90 -n 10 \
        -S -o graph/pggb/pggb_{wildcards.chromo}_{wildcards.grtype}'

        ln -s --relative graph/pggb/pggb_{wildcards.chromo}_{wildcards.grtype}/*.smooth.gfa {output}

        """





# rule cactus_combine_graph:


rule cactus_to_vg:
    input:"graph/cactus/cactus_{chromo}.hal"
    output:"graph/cactus/cactus_{chromo}.vg"
    threads: 20
    resources:
        mem_mb=2000,
        walltime="04:00"
    shell:
        """
        
        source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
        export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
        export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH

        hal2vg --inMemory {input} > {output} 

        """

rule cactus_to_gfa:
    input:"graph/cactus/cactus_{chromo}.vg"
    output:"graph/cactus/cactus_{chromo}.gfa"
    threads: 20
    resources:
        mem_mb=2000,
        walltime="04:00"
    shell:
        """
        vg convert -t {threads} -f {input} > {output}

        """

rule cactus_simplify:
    input:"graph/cactus/cactus_{chromo}.gfa"
    output:"graph/cactus/cactus_simple_{chromo}.gfa"
    threads: 20
    resources:
        mem_mb=2000,
        walltime="04:00"
    shell:
        """
        
        awk '$1 !~ /P/{{print;next}} $2 !~ /Anc/{{print}}' {input} > {output}


        """

rule cactus_stat:
    input:"graph/cactus/cactus_simple_{chromo}.gfa"
    output:"graph/cactus/cactus_stat_{chromo}.tsv"
    threads: 20
    resources:
        mem_mb=2000,
        walltime="04:00"
    params:
        ref=ref_genome
    shell:
        """
        
        ./graph_stat_mod.py -g {input} -o {output} -r {params.ref}.{wildcards.chromo}_{params.ref}

        """

localrules: cactus_combine_stat
rule cactus_combine_stat:
    input:expand("graph/cactus/cactus_stat_{chromo}.tsv",chromo=chromo_list),
    output:"graph/cactus/cactus_combine_stat.tsv"
    threads: 20
    resources:
        mem_mb=2000,
        walltime="04:00"
    run:

        from collections import defaultdict
        
        comb_stat=defaultdict(int)
        list_stat=defaultdict(list)

        outfile=open(output[0],"w")
        print(input)

        for statfile in input:
            with open(statfile) as infile:
                for line in infile:
                    print(line)
                    statid, statval = line.strip().split()
                    statval = int(statval)
                    comb_stat[statid] += statval
                    list_stat[statid].append(statval)

        for key,value in list_stat.items():
            value.append(comb_stat[key])
            print(key,*value,file=outfile)

        outfile.close()

#stat pggb and minigraph
# will delete later

rule pggb_stat:
    input:"graph/pggb_{chromo}/pggb_{chromo}.gfa"
    output:"graph/pggb/pggb_stat_{chromo}.tsv"
    threads: 20
    resources:
        mem_mb=2000,
        walltime="04:00"
    params:
        ref=ref_genome
    shell:
        """

        ./graph_stat_mod.py -g {input} -o {output} -r {wildcards.chromo}_{params.ref}


        """


localrules: pggb_combine_stat
rule pggb_combine_stat:
    input:expand("graph/pggb/pggb_stat_{chromo}.tsv",chromo=chromo_list),
    output:"graph/pggb/pggb_combine_stat.tsv"
    threads: 20
    resources:
        mem_mb=2000,
        walltime="04:00"
    run:

        from collections import defaultdict
        
        comb_stat=defaultdict(int)
        list_stat=defaultdict(list)

        outfile=open(output[0],"w")
        print(input)

        for statfile in input:
            with open(statfile) as infile:
                for line in infile:
                    print(line)
                    statid, statval = line.strip().split()
                    statval = int(statval)
                    comb_stat[statid] += statval
                    list_stat[statid].append(statval)

        for key,value in list_stat.items():
            value.append(comb_stat[key])
            print(key,*value,file=outfile)

        outfile.close()


rule modify_minigraph:
    input:  "graph/minigraph/graph_{chromo}.gfa"
    output: "graph/minigraph/graph_{chromo}_path.gfa"
    threads:10
    resources:
        mem_mb=2000 ,
        walltime= "04:00"
    shell:
        """
        vg convert -r 0 -g {input} -f > {output}

        """

rule minigraph_graph_statistics:
        input:"graph/minigraph/graph_{chromo}_path.gfa"
        output:"graph/minigraph/graph_{chromo}_stat.tsv"
        threads:10
        resources:
           mem_mb=8000 ,
           walltime= "04:00"
        params:
            ref=ref_genome
        shell:
           """

            ./graph_stat.py  -g {input}  -o {output} -r {wildcards.chromo}_{params.ref} -t minigraph

           """


localrules: minigraph_combine_stat
rule minigraph_combine_stat:
    input:expand("graph/minigraph/graph_{chromo}_stat.tsv",chromo=chromo_list),
    output:"graph/minigraph/minigraph_combine_stat.tsv"
    threads: 20
    resources:
        mem_mb=2000,
        walltime="04:00"
    run:

        from collections import defaultdict
        
        comb_stat=defaultdict(int)
        list_stat=defaultdict(list)

        outfile=open(output[0],"w")
        print(input)

        for statfile in input:
            with open(statfile) as infile:
                for line in infile:
                    print(line)
                    statid, statval = line.strip().split()
                    statval = int(statval)
                    comb_stat[statid] += statval
                    list_stat[statid].append(statval)

        for key,value in list_stat.items():
            value.append(comb_stat[key])
            print(key,*value,file=outfile)

        outfile.close()

