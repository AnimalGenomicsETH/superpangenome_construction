#!/usr/bin/env python

dirwork="/cluster/work/pausch/danang/psd/scratch/real_yak2"
chromo_list = list(range(25, 26))
assemb_list = ["UCD","Angus","Highland"]
ref_genome = "UCD"

rule all:
    input: expand("graph/minigraph/{chromo}/{chromo}_comb_coverage.tsv",chromo=chromo_list)

rule remap_graph:
    input:
        graph="graph/minigraph/graph_{chromo}.gfa",
        assembly="assembly/{chromo}/{breed}_{chromo}.fa"
    output:
        gfa=temp("graph/minigraph/{chromo}/{breed}_{chromo}_remap.gfa")
    threads: 10
    resources:
        mem_mb= 2000 ,
        walltime= "01:00"
    shell:
        """

        minigraph -t 10 --cov -x asm {input.graph} {input.assembly} > {output}

        """

rule process_coverage:
    input:rules.remap_graph.output
    output:"graph/minigraph/{chromo}/{breed}_{chromo}_coverage.tsv"
    threads: 10
    resources:
        mem_mb= 1000 ,
        walltime= "01:00"
    run:
        with open(input[0]) as infile, open(output[0],"w") as outfile:
            for line in infile:
                if line.startswith("S"):
                    token=line.strip().split()
                    node_id=token[1]
                    node_len=token[3].split(":")[2]
                    node_cov=token[-1].split(":")[2]
                    print(node_id,node_len,node_cov,file=outfile)

rule combine_coverage:
    input:expand("graph/minigraph/{{chromo}}/{breed}_{{chromo}}_coverage.tsv",breed=assemb_list)
    output:"graph/minigraph/{chromo}/{chromo}_comb_coverage.tsv"
    threads: 10
    resources:
        mem_mb= 1000 ,
        walltime= "01:00"
    run:
        import pandas as pd 
        
        
        for breed in assemb_list:
            infile=open(f"graph/minigraph/{wildcards.chromo}/{breed}_{wildcards.chromo}_coverage.tsv")
            if breed == ref_genome:
                combcov=pd.DataFrame([line.strip().split() for line in infile],
                                    columns=["node_id","node_len",breed])
            else:
                addcov=pd.DataFrame([line.strip().split() for line in infile],
                                    columns=["node_id","node_len",breed])
                
                combcov = pd.merge(combcov, addcov, on=["node_id","node_len"], how="outer")

            infile.close()
            
        covsimp=combcov.drop(columns=['node_id','node_len'])
        
        breed_list=assemb_list
        
        #determine color based on coverage
            
        def color_based_coverage(coverage,breed_list):
            color_list=[]

            for cov,breed in zip(coverage,breed_list):
                if float(cov) > 0:
                    color_list.append(breed)
            return ",".join(color_list)

        combcov['color'] = covsimp.apply(color_based_coverage,breed_list=breed_list,axis='columns')

        with open(output[0],"w") as outfile:
            combcov.to_csv(outfile,sep='\t',index=False)



        












