
##normalise direction of contigs
#rule download_fasta:
#    output:
#        expand('raw_assemblies/{sample}.fasta',sample=config['pangenome_samples'])
#    shell:
#        '''
#        ./prep_assemblies.sh
#        '''

rule split_fasta:
    input:
        'raw_assemblies/{sample}.fasta'
    output:
        fa = 'assemblies/{chromosome}/{sample}.fa',
        fai = 'assemblies/{chromosome}/{sample}.fa.fai',
    shell:
        '''
        echo "{wildcards.chromosome}" | seqtk subseq {input} - |\
        sed 's/>.*/>{wildcards.sample}/' > {output.fa}
        samtools faidx {output.fa}
        '''

rule repeatmasker_soft:
    input:
        #'assemblies/{chromosome}/{sample}.fa'
        '{path}/{sample}.fa'
    output:
        masked = temp('{path}/{sample}.fa.masked'),
        out = temp('{path}/{sample}.fa.out'),
        cat = temp('{path}/{sample}.fa.cat.gz'),
        tbl = temp('{path}/{sample}.fa.tbl'),
        #masked = temp('assemblies/{chromosome}/{sample}.fa.masked'),
        #out = temp('assemblies/{chromosome}/{sample}.fa.out'),
        #cat = temp('assemblies/{chromosome}/{sample}.fa.cat.gz')
    threads: 8
    resources:
        mem_mb = 1000,
        walltime = '24h'
    shell:
        '''
        RepeatMasker -pa $(({threads}/2)) -no_is -qq -xsmall \
        -lib {config[repeat_library]} {input}
        '''

rule mash_dist:
    input: 
        expand('assemblies/{chromosome}/{sample}.fa',chromosome=range(1,30),sample=config['pangenome_samples'])
    output:
        'communities'
    shell:
        '''
        cat {input} | bgzip -@ {threads} -c > $TMPDIR/genomes.fa.gz
        samtools faidx $TMPDIR/genomes.fa.gz
        mash dist -p 12 -s 5000 -i $TMPDIR/genomes.fa.gz $TMPDIR/genomes.fa.gz > bovinae.distances.tsv
        python paf2net.py -m bovinae.distances.tsv
        python net2communities.py -e bovinae.distances.tsv.edges.list.txt -w bovinae.distances.tsv.edges.weights.txt -n bovinae.distances.tsv.vertices.id2name.txt
        for i in {{0..29}}
          do
          chromosomes=$(cat bovinae.distances.tsv.edges.weights.txt.community.$i.txt | cut -f 3 -d '#' | sort | uniq | tr '\n' ' '
          echo "community $i --> $chromosomes"
        done
        '''

rule mash_triangle:
    input:
        lambda wildcards: expand('raw_assemblies/{sample}.fasta',sample=config[wildcards.sample_set])
    output:
        'tree/{sample_set}.lower_triangle.txt'
    threads: 6
    resources:
        mem_mb = 2500,
        walltime = '4h'
    shell:
        '''
        mash triangle -s 10000 -k 25 -p {threads} {input} | awk 'NR>1' > {output}
        '''

import numpy as np
from scipy.spatial.distance import squareform

def read_mash_triangle(mash_triangle,strip_leading=False):
    names, vals = [], []
    with open(mash_triangle,'r') as fin:
        for i,line in enumerate(fin):
            parts = line.rstrip().split()
            if strip_leading:
                names.append(parts[0].split('/')[1].split('.')[0])
            else:
                names.append(parts[0])
            vals.append(parts[1:]+[0])
    Q = np.asarray([np.pad(a, (0, len(vals) - len(a)), 'constant', constant_values=0) for a in vals],dtype=float)
    return names, (Q+Q.T)

def make_minigraph_order(mash_triangle,sequences=None):
    ref_ID = get_reference_ID()
    names, dists = read_mash_triangle(mash_triangle)
    sequences = sequences or names
    ref_ID_index = [i for i,s in enumerate(sequences) if ref_ID in s][0]

    return ' '.join(map(lambda k:k[0], sorted(zip(sequences,dists[ref_ID_index]),key=lambda k: k[1])))

rule minigraph_construct:
    input:
        mash_distances = expand(rules.mash_triangle.output,sample_set='pangenome_samples'),
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa', sample=pangenome_samples)
    output: 
        temp('graphs/minigraph/{chromosome}.basic.gfa')
    threads: 1
    resources:
        mem_mb = 20000,
        walltime = '24h'
    params:
        sample_order = lambda wildcards, input:make_minigraph_order(input.mash_distances[0],input.assemblies),
        L = config['minigraph']['L'],
        j = config['minigraph']['divergence']
    shell:
        '''
        minigraph -t {threads} -cxggs -j {params.j} -L {params.L} {params.sample_order} > {output}
        '''

rule minigraph_call:
    input:
        gfa = rules.minigraph_construct.output,
        sample = 'assemblies/{chromosome}/{sample}.fa'
    output:
        'graphs/minigraph/{chromosome}.{sample}.bed'
    params:
        L = config['minigraph']['L'],
        j = config['minigraph']['divergence']
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        '''
        minigraph -t {threads} -cxasm --call -j {params.j} -L {params.L} {input.gfa} {input.sample} > {output}
        '''

localrules: minigraph_path
rule minigraph_path:
    input:
        paths = expand('graphs/minigraph/{{chromosome}}.{sample}.bed',sample=filter(lambda x: x != get_reference_ID(), pangenome_samples)),
        gfa = 'graphs/minigraph/{chromosome}.basic.gfa'
    output:
        'graphs/minigraph/{chromosome}.gfa'
    params:
        samples = '\\n'.join(filter(lambda x: x != get_reference_ID(), pangenome_samples))
    shell:
        '''
        #needs special branch of mgutils.js
        {{ vg convert -r 0 -g {input.gfa} -f ; paste {input.paths} | mgutils.js path <(echo -e "{params.samples}") - | sed 's/s//g' ; }} > {output}
        '''

rule pggb_construct:
    input:
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa', sample=pangenome_samples),
    output:
        gfa = 'graphs/pggb/{chromosome}.gfa'
    threads: 16
    resources:
        mem_mb = 4000,
        walltime = '24h',
        scratch = '30G'
    params:
        pggb = config['pggb']['container'],
        _dir = lambda wildcards, output: Path(output[0]).with_suffix('').resolve(),
        fasta = '$TMPDIR/all.fa.gz',
        divergence = config['pggb']['divergence'],
        n_haplotypes = lambda wildcards, input: len(input.assemblies),
        segment_length = config['pggb']['segment_length']
    shell:
        '''
        cat {input} | bgzip -@ {threads} -c > {params.fasta}
        samtools faidx {params.fasta}
        
        mkdir -p {params._dir}

        singularity exec -B $TMPDIR -B {params._dir} {params.pggb} \
        pggb -i {params.fasta} -t {threads} \
        -s {params.segment_length} -p {params.divergence} -n {params.n_haplotypes} \
        --skip-viz --temp-dir $TMPDIR \
        -o {params._dir}

        mv {params._dir}/*.smooth.final.gfa {output.gfa}
        '''

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    if node.is_leaf():
        return f'{leaf_names[node.id]}:{parent_dist - node.dist}{newick}'
    else:
        if len(newick) > 0:
            newick = f'):{parent_dist - node.dist}{newick}'
        else:
            newick = ');'
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=f',{newick}')
        newick = f'({newick}'
        return newick

from scipy.cluster import hierarchy
localrules: cactus_seqfile
rule cactus_seqfile:
    input:
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa.masked', sample=pangenome_samples),
        mash_distances = expand(rules.mash_triangle.output,sample_set='pangenome_samples')
    output:
        temp('graphs/cactus/{chromosome}.seqfile')
    threads: 1
    resources:
        mem_mb = 1000
    run:
        names, dists = read_mash_triangle(input.mash_distances[0],True)
        Z = hierarchy.linkage(squareform(dists),method='average',optimal_ordering=True)
        tree = hierarchy.to_tree(Z, False)
        with open(output[0],'w') as fout:
            fout.write(get_newick(tree, tree.dist, names)+'\n')
            fout.write('\n'.join([f'{N} {P}' for (N,P) in zip(names,input.assemblies)]))

rule cactus_construct:
    input:
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa.masked', sample=pangenome_samples),
        seqFile = rules.cactus_seqfile.output
    output:
        jobStore = temp(directory('graphs/cactus/{chromosome}')),
        hal = temp('graphs/cactus/{chromosome}.hal')
    threads: 18
    resources:
        mem_mb = 3000,
        walltime = '4h',
        scratch = '50G'
    params:
        _asmDir = lambda wildcards, input: Path(input.assemblies[0]).parent.parent.resolve(),
        _dir = lambda wildcards, output: Path(output.hal).with_suffix('').resolve(),
        cactus = config['cactus']['container']
    shell:
        '''
        mkdir -p {params._dir}
        singularity exec -B $TMPDIR -B {params._dir} -B {params._asmDir} {params.cactus} \
        cactus --maxLocalJobs {threads} \
        --logLevel CRITICAL --workDir $TMPDIR \
        {output.jobStore}/jobStore {input.seqFile} {output.hal}
        '''

rule cactus_convert:
    input:
        rules.cactus_construct.output.hal
    output:
        'graphs/cactus/{chromosome}.gfa'
    threads: 1
    resources:
        mem_mb = 35000,
        walltime = '4h'
    params:
        _dir = lambda wildcards, output: Path(output[0]).parent.resolve(),
        cactus = config['cactus']['container']
    shell:
        '''
        singularity exec -B {params._dir} {params.cactus} \
        hal2vg --hdf5InMemory {input} |\
        vg convert -t {threads} -W -f - |\
        awk -v OFS='\t' '$1!~/P/ {{print;next}} $2!~/Anc/ {{split($2,a,"#");$2=a[1];print}}' |\
        sed 's/\*/0M/g' > {output}
        '''

rule odgi_stats:
    input:
        'graphs/{pangenome}/{chromosome}.gfa'
    output:
        'graphs/{pangenome}/{chromosome}.stats.yaml'
    threads: 1
    resources:
        mem_mb = 7500,
        walltime = '30'
    shell:
        '''
        odgi stats -i {input} -f -S -y -s -p > {output}
        '''

#LINE TO GET REF SIZE
#grep -hA 5 "HER"  graphs/minigraph/*.stats.yaml  | awk '$1~/nucleotides:/ {c+= $2} END {print c}'

rule summarise_gfa_stats:
    input:
        expand(rules.odgi_stats.output,chromosome=range(1,30),allow_missing=True)
    output:
        'graphs/{pangenome}/stats.yaml'
    shell:
        '''
        awk '{{ if ($1~/(length|edges|steps|file_size_in_bytes)/)  {{ A[$1]+=$2 }} else {{ if ($1~/nodes/&&A["nodes:"]==0) {{ A[$1]+=$2 }} }} }}  END {{ for (key in A) {{ print key,A[key] }} }}' {input} > {output}
        '''

rule odgi_flatten:
    input:
        'graphs/{pangenome}/{chromosome}.gfa'
    output:
        fasta = 'graphs/{pangenome}/{chromosome}.fa',
        bed = 'graphs/{pangenome}/{chromosome}.bed'
    resources:
        mem_mb = 7500
    shell:
        '''
        odgi flatten -i {input} -f {output.fasta} -b {output.bed}
        '''

localrules: repeat_summary
rule repeat_summary:
    input:
        expand('graphs/{pangenome}/{chromosome}.fa.tbl',pangenome=('minigraph','pggb','cactus'),chromosome=range(1,30))
    output:
        'graphs/repeats.csv'
    params:
        samples = lambda wildcards: '{' +','.join(list(pangenome_samples.keys())) +'}'
    shell:
        '''
        echo -e "pangenome\\tchromosome\\tSatellite\\tLTR\\tSINE\\tMasked\\tLINE\\tTotal" > {output}
        for i in {input}
        do
          mawk '/total|masked|SINE|LINE|Retroviral|Satellites/ {{ A[$1]+=$3 }} END {{ var=FILENAME; split (var,a,"/"); split(a[3],b,"."); printf a[2]"\\t"b[1]"\\t"; for (key in A) {{printf A[key]"\\t" }} {{printf "\\n"}} }}' $i >> {output}
        done
        for j in {{1..29}}
        do
          mawk -v C=$j '/total|masked|SINE|LINE|Retroviral|Satellites/ {{ A[$1]+=$3 }} END {{ printf "ASM\\t"C"\\t"; for (key in A) {{printf A[key]"\\t" }} {{printf "\\n"}} }}' assemblies/${{j}}/{params.samples}.fa.tbl >> {output}
          mawk -v C=$j '/total|masked|SINE|LINE|Retroviral|Satellites/ {{ A[$1]+=$3 }} END {{ printf "REF\\t"C"\\t"; for (key in A) {{printf A[key]"\\t" }} {{printf "\\n"}} }}' assemblies/${{j}}/HER.fa.tbl >> {output}
        done
        '''
