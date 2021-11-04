#!/usr/bin/env python

dirwork="/cluster/work/pausch/danang/psd/scratch/real_yak2"
#chromo_list = list(range(25, 26))
#chromo_list = [13,19,5,23,6,9]
#chromo_list = list(range(1,30))
#chromo_list=[2,3] + list(range(5,30))
chromo_list=list(range(1,3))
assemb_list = ["UCD","Angus","Highland","OBV","Brahman","Yak","BSW","Pied","Gaur","Nellore"]
#assemb_list = ["UCD","Angus"]
ref_genome = "UCD"
rep_library="BosTau9_repeat_library.fasta"
#temporary files for cactus
jobstore="/cluster/scratch/cdanang/cactus_dump"

rule all:
    input:
        #expand("assembly/{chromo}/{chromo}_rep/{breed}_{chromo}.fa.masked",chromo=chromo_list,breed=assemb_list),
        #expand("graph/cactus/cactus_{chromo}_seqfile.tsv",chromo=chromo_list),
        #expand("graph/cactus/cactus_{chromo}.hal",chromo=chromo_list),
	#expand("graph/cactus/cactus_{chromo}.vg",chromo=chromo_list),
	#expand("graph/cactus/cactus_{chromo}.gfa",chromo=chromo_list),
        #expand("graph/cactus/cactus_simple_{chromo}.gfa",chromo=chromo_list),
        #expand("graph/cactus/cactus_stat_{chromo}.tsv",chromo=chromo_list),
        #"graph/cactus/cactus_combine_stat.tsv",
        #"graph/pggb/pggb_combine_stat.tsv",
	#"graph/minigraph/minigraph_combine_stat.tsv",
        "graph/minigraph/graph_minigraph_combined.gfa",
        "graph/cactus/cactus_combined.gfa",
        expand("graph/cactus/cactus_{chromo}_snarl.tsv",chromo=chromo_list),
        expand("graph/minigraph/minigraph_{chromo}_snarl.tsv",chromo=chromo_list),


    
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
    input: expand("tree/{breed}.fa.msh", breed=assemb_list)
    output: "tree/combined_sketch.msh"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        mash paste {output} {input}

        """

rule estimate_distance:
    input: "tree/combined_sketch.msh"
    output: "tree/combined_distance.tsv"
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
        assemb=expand("assembly/{{chromo}}/{{chromo}}_rep/{breed}_{{chromo}}.fa.masked",breed=assemb_list),
        distance_file=rules.estimate_distance.output
    output:"graph/cactus/cactus_{chromo}_seqfile.tsv"
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
    input:"graph/cactus/cactus_{chromo}_seqfile.tsv"
    output:"graph/cactus/cactus_{chromo}.hal"
    threads: 20
    resources:
        mem_mb=5000,
        walltime="24:00"
    params:
        jobstore=jobstore + "/dump_{chromo}"
    shell:
        """


        source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
        export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
        export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH

        cactus --maxLocalJobs {threads} \
        {params.jobstore}  \
        {input} {output}

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



rule combine_cactus:
    input:expand("graph/cactus/cactus_simple_{chromo}.gfa",chromo=chromo_list)
    output:
        full_graph="graph/cactus/cactus_combined.vg",
        full_graph_gfa="graph/cactus/cactus_combined.gfa"
    threads: 20
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """

        vg combine {input} > {output.full_graph}
        vg convert -t {threads} -f {output.full_graph} > {output.full_graph_gfa}

        """

rule combine_minigraph:
    input:  expand("graph/minigraph/graph_{chromo}_path.vg",chromo=chromo_list)
    output:
        full_graph="graph/minigraph/graph_minigraph.vg",
        full_graph_gfa="graph/minigraph/graph_minigraph_combined.gfa",
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        vg combine {input} > {output.full_graph}

        vg convert -t {threads} -f {output.full_graph} > {output.full_graph_gfa}

        """



## deconstruct to identify snarl 

rule deconstruct_minigraph:
    input:"graph/minigraph/graph_{chromo}_path.gfa"
    output:"graph/minigraph/minigraph_{chromo}_snarl.tsv"
    threads: 10
    resources:
        mem_mb=2000,
        walltime="04:00"
    params:
        ref=ref_genome
    shell:
        """
        
        vg deconstruct -t {threads} -p {wildcards.chromo}_{params.ref} {input} > {output}
        
        """

rule deconstruct_cactus:
    input:"graph/cactus/cactus_simple_{chromo}.gfa"
    output:"graph/cactus/cactus_{chromo}_snarl.tsv"
    threads: 10
    resources:
        mem_mb=2000,
        walltime="04:00"
    params:
        ref=ref_genome
    shell:
        """
        # all level of snarls (including nested), -a
        # that only part of the paths, -e

        vg deconstruct -e -a -t {threads} -p {params.ref}.{wildcards.chromo}_{params.ref} {input} > {output}

        """


rule deconstruct_pggb:
    input:"graph/cactus/cactus_simple_{chromo}.gfa"
    output:"graph/cactus/cactus_{chromo}_snarl.tsv"
    threads: 10
    resources:
        mem_mb=2000,
        walltime="04:00"
    params:
        ref=ref_genome
    shell:
        """
        # all level of snarls (including nested), -a
        # that only part of the paths, -e

        vg deconstruct -e -a -t {threads} -p {wildcards.chromo}_{params.ref} {input} > {output}

        """
