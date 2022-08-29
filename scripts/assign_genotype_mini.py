#!/usr/bin/env python
import re
from collections import defaultdict 

assemb_dir="/cluster/work/pausch/danang/psd/scratch/assembly"
breed_list="BSW,Pied,Highland,Angus,OBV,Simmental,Brahman,Nellore,Gaur,Bison,Yak"
breed_list=breed_list.split(",")
chromo_list=range(1,30)

rule all:
    input:
         # expand("snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.bed",breed=breed_list,chromo=chromo_list),
         # expand("snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}_path.gfa",breed=breed_list,chromo=chromo_list),
         # expand("snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.vcf",breed=breed_list,chromo=chromo_list),
         expand("snarl/snarl_new/minigraph/call_sv/comb/comb_{chromo}.gfa",chromo=chromo_list),
         expand("snarl/snarl_new/minigraph/call_sv/comb/comb_{chromo}.vcf",chromo=chromo_list)



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
                        P_line[path_id].extend(path_node)
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
                path_join=",".join(value)
                print(f"P\t{wildcards.chromo}_{wildcards.breed}_{key}\t{path_join}\t*\t",file=outfile)
        


rule convert_vg:
    input:
         "snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.gfa"
    output:
         "snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}_path.gfa"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        vg convert -r 0 -g {input} -f > {output}

        """


rule deconstruct:
    input:
        "snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}_path.gfa"
    output:
        "snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{chromo}.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    shell:
        """

        vg deconstruct -p {wildcards.chromo}_UCD -d 1 -a -e {input} > {output}

        """

rule add_path_combined:
    input:
        call=expand("snarl/snarl_new/minigraph/call_sv/{breed}/{breed}_{{chromo}}.bed",breed=breed_list),
        graph="graph/minigraph/minigraph_{chromo}.gfa"
    output:
        "snarl/snarl_new/minigraph/call_sv/comb/comb_{chromo}.gfa"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    run:
         
        P_comb=[]
        


        for input_file in input.call:
            breed=input_file.split("/")[-1].split("_")[0]
            chromo=input_file.split("/")[-1].split("_")[1].split(".")[0]
            P_line=defaultdict(list)
            path_id=1

            with open(input_file) as infile:
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
                        P_line[f"{breed}_{chromo}_{path_id}"].append(start_node)
                        if path_node:
                            P_line[f"{breed}_{chromo}_{path_id}"].extend(path_node)
                        P_line[f"{breed}_{chromo}_{path_id}"].append(end_node)
                        prev_node=end_node
                    else:
                        if prev_node==start_node:
                            if path_node:
                                P_line[f"{breed}_{chromo}_{path_id}"].extend(path_node)
                            P_line[f"{breed}_{chromo}_{path_id}"].append(end_node)
                        else:
                            path_id += 1
                            P_line[f"{breed}_{chromo}_{path_id}"].append(start_node)
                            if path_node:
                                P_line[f"{breed}_{chromo}_{path_id}"].extend(path_node)
                            P_line[f"{breed}_{chromo}_{path_id}"].append(end_node)
                        prev_node=end_node
                P_comb.append(P_line)


        with open(input.graph) as infile, open(output[0],"w") as outfile:
            for line in infile:
                print(line.strip(),file=outfile)
            for comp in P_comb:
                for key,value in comp.items():
                    path_join=",".join(value)
                    print(f"P\t{key}\t{path_join}\t*\t",file=outfile)
        
rule convert_vg_comb:
    input:
        "snarl/snarl_new/minigraph/call_sv/comb/{chromo}_comb.gfa"
    output:
        "snarl/snarl_new/minigraph/call_sv/comb/{chromo}_comb_path.gfa"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        vg convert -r 0 -g {input} -f > {output}

        """


rule deconstruct_comb:
    input:
        "snarl/snarl_new/minigraph/call_sv/comb/{chromo}_comb_path.gfa"
    output:
        "snarl/snarl_new/minigraph/call_sv/comb/{chromo}_comb.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "01:00"
    shell:
        """

        vg deconstruct -p {wildcards.chromo}_UCD -d 1 -a -e {input} > {output}

        """
