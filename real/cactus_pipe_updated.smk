#!/usr/bin/env python
scrdir="/cluster/scratch/cdanang/scr_cac"
#make it chromosome-by-chromosome 
#this use input GFA from minigraph, aligned, and masked in a single pre-process step 


localrules: create_seq_file 
rule create_seq_file:
        input:
            graph="graph/minigraph/graph_{chromo}.gfa",
            fasta=expand("assembly/{{chromo}}/{anim}_{{chromo}}.fa",anim=assemb_list)
        output:
            seqfile="graph/cactus/{chromo}/seqfile_{chromo}.txt",
            seqfile_brnn="graph/cactus/{chromo}/seqfile_{chromo}_brnn.txt",
            seqfile_masked="graph/cactus/{chromo}/seqfile_{chromo}_secmask.txt"
        params:
            prefix="graph/cactus/{chromo}",
            anims=assemb_list,
            chromo="{chromo}"
        shell:
           """

           for anim in {params.anims}
           do
                echo {params.chromo}_${{anim}} assembly/{params.chromo}/${{anim}}_{params.chromo}.fa >> {params.prefix}/seqfile_{params.chromo}.txt
                echo {params.chromo}_${{anim}} {params.prefix}/${{anim}}_{params.chromo}_brnn.fa >> {params.prefix}/seqfile_{params.chromo}_brnn.txt
                echo {params.chromo}_${{anim}} {params.prefix}/${{anim}}_{params.chromo}_secmask.fa >> {params.prefix}/seqfile_{params.chromo}_secmask.txt
           done 

           """

#preprocess with dna-brnn
#map with graphmap
#preprocess with paf from graphmap to add in gaps

#Step 1: Masked with DNA_BRNN
rule cactus_preprocess_brnn:
        input:
            seqfile="graph/cactus/{chromo}/seqfile_{chromo}.txt",
            seqfile_masked="graph/cactus/{chromo}/seqfile_{chromo}_brnn.txt"
        output:expand("graph/cactus/{{chromo}}/{anim}_{{chromo}}_brnn.fa",anim=assemb_list)
        threads:32 
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        params:
            chromo = "{chromo}"
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH


           cactus-preprocess $PWD/{params.chromo}_src_masked \
           {input.seqfile} {input.seqfile_masked} \
           --realTimeLogging --binariesMode local \
           --maskAlpha --brnnCores 8


           """

# Step 2: Masking with coverage map of PAF
# Step 2.1: Map brnn_fasta to the graph

rule cactus_preprocess_graphmap:
        input: 
            graph="graph/minigraph/graph_{chromo}.gfa",
            fasta=expand("graph/cactus/{{chromo}}/{anim}_{{chromo}}_brnn.fa",anim=assemb_list),
            seqfile="graph/cactus/{chromo}/seqfile_{chromo}_brnn.txt",
        output: 
            paf="graph/cactus/{chromo}/cactus_{chromo}.paf"
        threads: 10
        resources:
           mem_mb= 5000,
           walltime= "04:00"
        params:
            prefix="graph/cactus/{chromo}",
            scrdir=scrdir +"/jobstore_map_{chromo}",
            anims=assemb_list,
            chromo="{chromo}"
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH


            cactus-graphmap $PWD/{params.chromo}_src  \
            {input.seqfile} \
            {input.graph} {output.paf} \
            --outputFasta {params.prefix}/cactus_{wildcards.chromo}_sec.fa \
            --realTimeLogging

           """

#Step 2.2: Use PAF coverage for final masking 
rule cactus_preprocess_second_masked:
        input:
            seqfile="graph/cactus/{chromo}/seqfile_{chromo}_brnn.txt",
            seqfile_masked="graph/cactus/{chromo}/seqfile_{chromo}_secmask.txt",
            paf="graph/cactus/{chromo}/cactus_{chromo}.paf"
        output:expand("graph/cactus/{{chromo}}/{anim}_{{chromo}}_secmask.fa",anim=assemb_list)
        threads:32 
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        params:
            chromo = "{chromo}"
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH


           cactus-preprocess $PWD/{params.chromo}_src_masked \
           {input.seqfile} {input.seqfile_masked} \
           --realTimeLogging --binariesMode local \
           --maskFile {input.paf} \
           --minLength 100000
           


           """


#Step 3: Real cactus alignment
rule cactus_align:
        input:
              seqfile_masked="graph/cactus/{chromo}/seqfile_{chromo}_brnn.txt",
              paf="graph/cactus/{chromo}/cactus_{chromo}.paf",
              fasta=expand("graph/cactus/{{chromo}}/{anim}_{{chromo}}_secmask.fa",anim=assemb_list)
        output:
              gfa="graph/cactus/{chromo}/cactus_{chromo}.gfa.gz",
              hal="graph/cactus/{chromo}/cactus_{chromo}.hal",
              vg="graph/cactus/{chromo}/cactus_{chromo}.vg"
        threads:32
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        params: 
           ref = ref_genome
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH

           cactus-align $PWD/{wildcards.chromo}_src_align \
           {input.seqfile_masked} \
           {input.paf} {output.hal} \
           --pangenome --pafInput \
           --realTimeLogging \
           --outGFA --outVG --reference {wildcards.chromo}_{params.ref}

           """

rule cactus_drop_paths:
        input:"graph/cactus/{chromo}/cactus_{chromo}.vg"
        output:"graph/cactus/{chromo}/cactus_drop_{chromo}.vg"
        threads:10
        resources:
           mem_mb=2000 ,
           walltime= "04:00"
        params:
        shell:
           """
            
          vg paths -Q _MINIGRAPH -d -v {input} > {output}

           """

rule cactus_combine:
        input:expand("graph/cactus/{chromo}/cactus_drop_{chromo}.vg",chromo=chromo_list)
        output:"graph/cactus/graph_cactus_combined.vg"
        threads:10
        resources:
           mem_mb= 5000,
           walltime= "04:00"
        params:
        shell:
           """
            
            vg combine {input} > {output}

           """

rule cactus_convert_gfa:
        input:"graph/cactus/graph_cactus_combined.vg"
        output:"graph/cactus/graph_cactus_combined.gfa"
        threads:10
        resources:
           mem_mb= 5000,
           walltime= "04:00"
        shell:
           """

           vg convert -t {threads} -f {input} > {output}
           """

rule cactus_graphmap_join:
        input:
            hal=expand("graph/cactus/{chromo}/cactus_{chromo}.hal",chromo=chromo_list),
            vg=expand("graph/cactus/{chromo}/cactus_{chromo}.vg",chromo=chromo_list)
        output:touch("combine_finished.tsv")
        threads:32
        resources:
           mem_mb=2000 ,
           walltime= "04:00"
        shell:
           """

           source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
           export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
           export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH

           cactus-graphmap-join $PWD/src_combine --outDir $PWD/graph/cactus_combine \
           --outName cactus_comb --reference UCD --realTimeLogging --clipLength 100000 --wlineSep . \
           --vg {input.vg} 
           """

rule chop_cactus:
        input:"graph/cactus/graph_cactus_combined.vg"
        output:"graph/cactus/graph_cactus_chop.vg"
        threads:10
        resources:
           mem_mb= 5000,
           walltime= "04:00"
        shell:
           """

            vg mod -t {threads} -X 1000 {input} > {output}

           """

rule convert_cactus_chop_gfa:
        input:"graph/cactus/graph_cactus_chop.vg"
        output:"graph/cactus/graph_cactus_chop.gfa"
        threads:10
        resources:
           mem_mb= 5000,
           walltime= "04:00"
        shell:
           """
            vg convert -f {input} > {output}
            
           """
