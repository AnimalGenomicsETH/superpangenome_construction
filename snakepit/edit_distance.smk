import subprocess
import tempfile

rule seperate_centro_telo_regions:
    input:
        repeats = 'assemblies/{chromosome}/{sample}.fa.out',
        fasta = 'assemblies/{chromosome}/{sample}.fa',
        fai = 'assemblies/{chromosome}/{sample}.fa.fai'
    output:
        trimmed = temp('edit_distance/{sample}.{chromosome}.chunks.trimmed.fasta'),
        untrimmed = temp('edit_distance/{sample}.{chromosome}.chunks.untrimmed.fasta')
    params:
        centromere_min_score = 50000,
        telomere_min_score = 1000,
        window = config['graphaligner_parameters']['window'],
        chr_size = lambda wildcards, input: int([l.split()[1] for l in open(input.fai)][0])
    resources:
        mem_mb = 500,
        walltime = '60',
        disk_scratch = 1
    run:
        centro_region = subprocess.run(f'awk \'/Satellite\/centr/ {{print $5"\\t"$6"\\t"$7"\\t"$1}}\' {input.repeats} | bedtools merge -d 250000 -i - -c 4 -o sum | head -n 1',shell=True,capture_output=True).stdout.decode("utf-8")
        telo_region = subprocess.run(f'awk \'/\\(TTAGGG\\)n/||/\\(TAGGGT\\)n/||/\\(AGGGTT\\)n/||/\\(GGGTTA\\)n/||/\\(GGTTAG\\)n/||/\\(GTTAGG\\)n/ {{print $5"\\t"$6"\\t"$7"\\t"$1}}\' {input.repeats} | bedtools merge -d 1000 -i - -c 4 -o sum | tail -n 1',shell=True,capture_output=True).stdout.decode("utf-8")

        start = 1 if (not centro_region or int(centro_region.split()[-1]) < params.centromere_min_score) else int(centro_region.split()[2])
        end = params.chr_size if (not telo_region or int(telo_region.split()[-1]) < params.telomere_min_score) else int(telo_region.split()[1])
        
        masked_regions = ''
        if start != 1:
            masked_regions += '\n'.join((f'{wildcards.sample}:{i}-{min(start,i+params.window-1)}' for i in range(1,start+1,params.window)))
        if end != params.chr_size:
            if masked_regions:
                masked_regions += '\n'
            masked_regions += '\n'.join((f'{wildcards.sample}:{i}-{min(params.chr_size,i+params.window-1)}' for i in range(end,params.chr_size+1,params.window)))

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            temp.write(masked_regions)
            temp.seek(0)
            subprocess.run(f'samtools faidx -r {temp.name} {input.fasta} > {output.untrimmed}',shell=True)
    
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            temp.write('\n'.join((f'{wildcards.sample}:{i}-{min(end,i+params.window-1)}' for i in range(start,end+1,params.window))))
            temp.seek(0)
            subprocess.run(f'samtools faidx -r {temp.name} {input.fasta} > {output.trimmed}',shell=True)

        #with open(output[0],'w') as fout:
        #    fout.write(f'{wildcards.sample},{wildcards.chromosome},{start-1},{end-start-1},{params.chr_size-end}')

def get_threads(wildcards,input):
    if input.size_mb == 0:
        return 1
    if wildcards.trimmed == 'trimmed':
        return 4
    else:
        return 1
    if pangenome == 'minigraph':
        return 2
    elif pangenome == 'cactus':
        return 2
    elif pangenome == 'pggb':
        return 2

def get_memory(wildcards,input):
    if input.size_mb == 0:
        return 100
    if wildcards.trimmed == 'trimmed':
        return 15000
    else:
        return 30000

rule graphaligner:
    input:
        gfa = 'graphs/{pangenome}/{chromosome}.gfa',
        fasta = lambda wildcards: 'edit_distance/{sample}.{chromosome}.chunks.{trimmed}.fasta'
    output:
        temp('edit_distance/{sample}.{chromosome}.{pangenome}.{preset}.{trimmed}.gaf')
    params:
        preset = lambda wildcards: config['graphaligner_parameters'][wildcards.preset]
    threads: 4#lambda wildcards,input: get_threads(wildcards,input)
    resources:
        mem_mb = lambda wildcards,input: get_memory(wildcards,input),
        walltime = '4:00'
    shell:
        '''
        if [ -s {input.fasta} ]; then
          GraphAligner -g {input.gfa} -f {input.fasta} --multimap-score-fraction 1 -a {output} -t {threads} --discard-cigar {params.preset}
        else
          touch {output}
        fi
        '''

localrules: bin_edit_distance
rule bin_edit_distance:
    input:
        rules.graphaligner.output[0]
    output:
        'edit_distance/{sample}.{chromosome}.{pangenome}.{preset}.{trimmed}.dist'
    shell:
        '''
         awk '{{M[$1]+=$10;L[$1]+=$11}} END {{ for (key in M) {{ print key,M[key],L[key] }} }}' {input} | sort -k1,1V > {output}
        '''

localrules: gather_edit
rule gather_edit:
    input:
        rules.bin_edit_distance.output[0]
    output:
        'edit_distance/{sample}.{chromosome}.{pangenome}.{preset,strict|lenient}.{trimmed}.stat'
    shell:
        '''
        awk '{{split($1,a,":");split(a[2],b,"-"); D+=(b[2]-b[1]);L+=$2;M+=$3}} END {{print L/M,M/D,L,M,D}}' {input} > {output}
        '''

#for p in minigraph cactus; do for i in {1..29}; do for j in Angus Bison Brahman BSW Gaur Highland Nellore OBV Pied Simmental UCD Yak; do for k in trimmed untrimmed; do echo $p $j $i $k $(awk '{{split($1,a,":");split(a[2],b,"-"); D+=(b[2]-b[1]);L+=$2;M+=$3}} END {{print L,M,D}}' edit_distance/${i}_${j}.${p}.${k}.dist); done;done;done;done > big.csv
