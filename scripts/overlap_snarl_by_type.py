#!/usr/bin/env python

#chromo_list=list(range(1,30))
chromo_list=range(1,2)
#prog_list=["cactus","pggb","minigraph"]
prog_list=["cactus","pggb"]
#chromo_list.remove(5)
#chromo_list.remove(6)
breed_list="BSW,OBV,Pied,Highland,Angus,UCD,Simmental,Brahman,Nellore,Gaur,Bison,Yak"
breed_list=breed_list.split(",")


rule all:
    input:
        # expand("{prog}/{prog}_{chromo}_sv.fa",chromo=chromo_list),
        # expand("{prog}/{prog}_{chromo}_sv.fa.tbl",chromo=chromo_list),
        # expand("{prog}/{prog}_{chromo}_sv.fa.tbl",chromo=chromo_list),
        # expand("{prog}/{prog}_{chromo}_sv_annot.vcf",chromo=chromo_list,prog=prog_list),
        # expand("{prog}/{prog}_{chromo}_sv_norep.vcf",chromo=chromo_list,prog=prog_list),
        # expand("merged/{chromo}_merged_norep.tsv",chromo=chromo_list),
        # expand("merged/{chromo}_merged_norep.vcf",chromo=chromo_list),
        # expand("merged/{chromo}_merged_all.vcf",chromo=chromo_list),
        expand("{prog}/{breed}/{prog}_{chromo}_{breed}_sv.vcf",prog=prog_list,chromo=chromo_list,breed=breed_list),
        expand("merged/{breed}/{chromo}_{breed}_merged.vcf",breed=breed_list,chromo=chromo_list)

rule trim_allele:
    input:
        "{prog}/{prog}_{chromo}_snarl.vcf"
    output:
        "{prog}/{prog}_{chromo}_sv_trim2.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        if [[ {wildcards.prog} == "minigraph" ]]; then 
        
            ln -r -s {input} {output}

        else 

            bcftools view -a {input} |  bcftools view -i 'abs(ILEN)>=50' > {output}
    
        fi 

        """


rule extract_repeat:
    input:
        "{prog}/{prog}_{chromo}_sv_trim2.vcf"
    output:
        "{prog}/{prog}_{chromo}_sv.fa"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"

    run: 
        with open(input[0]) as infile, open(output[0],"w") as outfile:
            for line in infile:
                if not line.startswith("#"):
                    token=line.strip().split()
                    varid=f"{token[0]}_{token[1]}"
                    all_allele=[token[3]]
                    all_allele.extend(token[4].split(","))
                    all_allele.sort(key=lambda x:len(x))
                    max_allele=all_allele[-1]
                    print(f">{varid}")
                    print(max_allele)

    # shell:
        # """
        # ./extract_fasta_len.py {input} > {output}
        # """

rule repeat_masker:
    input:
         "{prog}/{prog}_{chromo}_sv.fa"
    output:
         "{prog}/{prog}_{chromo}_sv.fa.out"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        RepeatMasker -qq -no_is -lib /cluster/work/pausch/alex/Libraries/BosTau9_repeat_library.fasta {input}

        """

rule annotate_var:
    input:
        vcf="{prog}/{prog}_{chromo}_sv_trim2.vcf",
        repfile="{prog}/{prog}_{chromo}_sv.fa.out"
    output:
        "{prog}/{prog}_{chromo}_sv_annot.vcf"
    threads: 10 
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    run:
        from collections import defaultdict
        
        repcol=defaultdict(set)
        outfile=open(output[0],"w")
        #get dict of the rep 
        with open(input.repfile) as infile:
            next(infile)
            next(infile)
            next(infile)

            for line in infile:
                token=line.strip().split()
                repcol[token[4]].add(token[10])

            with open(input.vcf) as infile:
                for line in infile:
                    if not line.startswith("#"):
                        token=line.strip().split()
                        varid=f"{token[0]}_{token[1]}"
                        repid=repcol.get(varid,["norep"])
                        token[6]=",".join(repid)
                        print(*token,sep="\t",file=outfile)


rule remove_norep:
    input:
        header="{prog}/{prog}_{chromo}_sv_trim2.vcf",
        annot="{prog}/{prog}_{chromo}_sv_annot.vcf"
    output:
        "{prog}/{prog}_{chromo}_sv_norep.vcf"
    threads: 2
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        cat <(grep "#" {input.header}) <(grep -v "#" {input.annot} | grep "norep" )  |
        sed 's/_UCD//g' > {output}

        """


localrules: create_jasmine_file
rule create_jasmine_file:
    input:
        expand("{prog}/{prog}_{{chromo}}_sv_norep.vcf",prog=prog_list)
    output:
        "merged/{chromo}_merged_norep.tsv"
    shell:
        """
        echo {input} | tr ' ' '\n' > {output}

        """


rule merged_jasmine:
    input:
        "merged/{chromo}_merged_norep.tsv"
    output:
        "merged/{chromo}_merged_norep.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        #conda activate jasmine 

        jasmine file_list={input} \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_read=0 \
        --ignore_strand --threads=8 --allow-intersample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output}

        """


rule overlap_all:
    input:
        expand("{prog}/{prog}_{{chromo}}_sv_annot.vcf",prog=prog_list)
    output:
        listfiles="merged/{chromo}_merged_all_list.tsv",
        merged="merged/{chromo}_merged_all.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        echo {input} | tr ' ' '\n' > {output.listfiles}
        

        jasmine file_list={output.listfiles} \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_reads=0 \
        --ignore_strand --threads=8 --allow_intrasample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output.merged}

        """


#rule overlap assemblies

rule persample_var:
    input:
        "{prog}/{prog}_{chromo}_sv_trim2.vcf"
    output:
        header="{prog}/header/{chromo}_{prog}_{breed}_newheader.tsv",
        corvcf="{prog}/{breed}/{prog}_{chromo}_svcor.vcf",
        breedsv="{prog}/{breed}/{prog}_{chromo}_{breed}_sv.vcf"
    threads: 5
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        awk '$1 ~/#/' {input} | 
        awk '$0 !~ /ID=AT/{{print $0;next}}{{printf "%s\\n##INFO=<ID=CONFLICT,Number=1,Type=String,Description=Uncertain Paths>\\n", $0}}' |
        sed 's/_UCD//g' > {output.header}
        
        bcftools reheader -h {output.header} {input} > {output.corvcf}

        bcftools view --min-ac 1 -I -a -s {wildcards.chromo}_{wildcards.breed} {output.corvcf} > {output.breedsv}

        """


rule overlap_breed:
    input:
        expand("{prog}/{{breed}}/{prog}_{{chromo}}_{{breed}}_sv.vcf",prog=prog_list)
    output:
        listfiles="merged/{breed}/{chromo}_{breed}_list.tsv",
        merged="merged/{breed}/{chromo}_{breed}_merged.vcf"
    threads: 10
    resources:
        mem_mb= 2000,
        walltime= "04:00"
    shell:
        """

        echo {input} | tr ' ' '\n' > {output.listfiles}
        

        jasmine file_list={output.listfiles} \
        genome_file=/cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa \
        spec_reads=0 \
        --ignore_strand --threads=8 --allow_intrasample --dup_to_ins --output_genotypes --normalize_type --pre_normalize \
        out_file={output.merged}

        """

