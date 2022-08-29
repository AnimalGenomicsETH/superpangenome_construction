#!/usr/bin/env python
import re
from collections import defaultdict 

assemb_dir="/cluster/work/pausch/danang/psd/scratch/assembly"
breed_list=["Angus"]
chromo_list=[25]

rule all:
    input:
         expand("snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.bed",breed=breed_list,chromo=chromo_list),
         expand("snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.gfa",breed=breed_list,chromo=chromo_list)


rule realign:
    input:
        graph="graph/minigraph/minigraph_{chromo}.gfa",
        assemb= assemb_dir + "/{chromo}/{breed}_{chromo}.fa"
    output:
        "snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.bed"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    shell:
        """

        minigraph19 -cxasm --call {input.graph} {input.assemb} > {output}

        """


# rule match_geno:
    # input:
        # snarl="snarl/snarl_new/minigraph/minigraph_{chromo}_snarl.vcf",
        # call="snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.bed"
    # output:
        # "snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}_match.vcf"
    # threads: 10
    # resources:
        # mem_mb= 2000,
        # walltime= "04:00"
    # run:

        # callgeno=dict()
        # with open(input.call) as infile:
            # for line in infile:
                # token=line.strip().split()
                # varid=f"{token[3]}_{token[4]}"
                # geno=token[-1].split(":")[0]
                # callgeno[varid]=geno
        

        # with open(input.snarl) as infile,open():
            # for line in infile:
                # if line



rule add_path:
    input:
        call="snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.bed",
        graph="graph/minigraph/minigraph_{chromo}.gfa"
    output:
        "snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.gfa"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    run:
         
        P_line=defaultdict(list)
        path_id=1

        with open(input.call) as infile:
            for ind,line in enumerate(infile):
                token=line.strip().split()
                start=token[3]
                start_node=f"{start[1:]}+" if start[0] == ">" else f"{start[1:]}-"
                end=token[4]
                end_node=f"{end[1:]}+" if end[0] == ">" else f"{end[1:]}-"
                path=token[-1].split(":")[0]
                path_node=[]
                if path == "." and path == "*":
                    path_node=""
                else:
                    all_nodes=re.findall("[>|<]s[0-9]+",path)
                    for node in all_nodes:  
                        path_node.append(f"{node[1:]}+" if node[0] == ">" else f"{node[1:]}-")

                if not ind:
                    P_line[path_id].append(start_node)
                    if path_node:
                        P_line.extend(path_node)
                    P_line[path_id].append(end_node)
                    prev_node=end_node
                else:
                    if prev_node==start_node:
                        if path_node:
                            P_line[path_id].extend(path_node)
                        P_line[path_id].append(end_node)
                    else:
                        path_id += 1
                        P_line[path_id].append(start_node)
                        if path_node:
                            P_line[path_id].extend(path_node)
                        P_line[path_id].append(end_node)
                    prev_node=end_node


        with open(input.graph) as infile, open(output[0],"w") as outfile:
            for line in infile:
                print(line.strip(),file=outfile)
            for key,value in P_line.items():
                path_join=",".join(value).replace("s","")
                print(f"P\t25_Angus_{key}\t{path_join}\t",file=outfile)
        


                
