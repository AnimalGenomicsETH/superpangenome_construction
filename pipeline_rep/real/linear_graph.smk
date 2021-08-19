## linear mapping 

breed_list=["UCD","OBV"]
anim_list=["BSWCHEM110294048847"]
fastq_dir="/cluster/work/pausch/inputs/fastq/BTA"


rule all:
    input: expand("mapped/linear_{breeds}_{anim}_mapping_stat.tsv",breeds=breed_list,anim=anim_list)

rule construct_linear:
        input:"assembly/{breeds}_aut.fa"
        output:
            vg="graph/linear/{breeds}_linear.vg",
            gfa="graph/linear/{breeds}_linear.gfa",
            xg="graph/linear/{breeds}_linear.xg"
        threads:32
        resources:
           mem_mb=2000 ,
           walltime= "01:00"
        shell:
           """
           vg construct -m 1000 -t {threads} -r {input} > {output.vg}

           vg convert -t {threads} -f {output.vg} > {output.gfa} 

           vg convert -t {threads} -x {output.vg} > {output.xg}

           """

rule create_gbwt_linear:
        input: "graph/linear/{breeds}_linear.xg"
        output:  "graph/linear/{breeds}_linear.giraffe.gbwt"
        threads: 32
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        shell:
           """

           vg gbwt --index-paths -x {input} -o {output}

           """

rule identify_snarl:
        input: "graph/linear/{breeds}_linear.xg"
        output: "graph/linear/{breeds}_linear.snarl"
        threads: 32
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        shell:
           """

           vg snarls -t {threads} -T {input} > {output}

           """


rule create_distance_linear:
        input: xg="graph/linear/{breeds}_linear.xg",
               snarl="graph/linear/{breeds}_linear.snarl"
        output:"graph/linear/{breeds}_linear.dist"
        threads: 32
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        shell:
           """

           vg index -t {threads} {input.xg} \
            -s {input.snarl} -j {output}

           """

rule create_gbwt_graph_linear:
        input: 
            xg="graph/linear/{breeds}_linear.xg",
            gbwt="graph/linear/{breeds}_linear.giraffe.gbwt"
        output: "graph/linear/{breeds}_linear.gg"
        threads: 32
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        shell:
           """

           vg gbwt -g {output} -x {input.xg} {input.gbwt}

           """

rule create_index_minimizer:
        input:
            xg="graph/linear/{breeds}_linear.xg",
            gbwt="graph/linear/{breeds}_linear.giraffe.gbwt"
        output: "graph/linear/{breeds}_linear.min"
        threads: 32
        resources:
           mem_mb= 2000 ,
           walltime= "04:00"
        shell:
           """

           vg minimizer -t {threads} -g {input.gbwt} -i {output} {input.xg}

           """
rule giraffe_mapping:
        input:
            xg="graph/linear/{breeds}_linear.xg",
            gg="graph/linear/{breeds}_linear.gg",
            gbwt="graph/linear/{breeds}_linear.giraffe.gbwt",
            mini="graph/linear/{breeds}_linear.min",
            dist="graph/linear/{breeds}_linear.dist",
            r1=fastq_dir + "/{anim}_R1.fastq.gz",
            r2=fastq_dir + "/{anim}_R2.fastq.gz"
        output:"mapped/linear_{breeds}_{anim}.gam"
        threads:32
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        params:
            fastq_dir=fastq_dir
        shell:
           """

           vg giraffe -t {threads} \
           -x {input.xg} \
           -g {input.gg} \
           -H {input.gbwt} \
           -m {input.mini} \
           -d {input.dist} \
           -p -f {input.r1} -f {input.r2} > {output}

           """

rule graph_mapping_statistics:
        input:"mapped/linear_{breeds}_{anim}.gam"
        output:"mapped/linear_{breeds}_{anim}_mapping_stat.tsv"
        threads:10
        resources:
           mem_mb= 2000,
           walltime= "04:00"
        shell:
           """
        
            vg stats -a {input} > {output}
           
           """   