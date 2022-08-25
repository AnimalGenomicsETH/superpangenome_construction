import subprocess
import tempfile
from pathlib import PurePath

breeds = ['Angus', 'Bison', 'Brahman', 'BSW', 'Gaur', 'Highland', 'Nellore', 'OBV', 'Pied', 'Simmental', 'UCD', 'Yak', 'BS1','BS2','BS3','BS4','BS5','BS6','BS7','BS8']
HiFi_breeds = ['BSW','Gaur','Nellore','OBV','Pied','BS1','BS2','BS3','BS4','BS5','BS6','BS7','BS8']

localrules: bin_edit_distance

rule all:
    input:
        expand('edit_distance2/{chrom}_{asm}.{pangenome}.untrimmed.dist',chrom=range(1,30),asm=breeds,pangenome=('pggb','cactus','minigraph')),
        expand('edit_distance2/{chrom}_{asm}.{pangenome}.trimmed.dist',chrom=range(1,30),asm=breeds,pangenome=('pggb','cactus','minigraph'))


rule RepeatMasker:
    input:
        '/cluster/work/pausch/danang/psd/scratch/real_comp/graph/map_lr/splitfa/{chrom}/{chrom}_BS{N}.fa'
    output:
        '/cluster/work/pausch/danang/psd/scratch/assembly/{chrom}/{chrom}_rep/{chrom}_BS{N}.fa.out'
    params:
        library = '/cluster/work/pausch/alex/Libraries/BosTau9_repeat_library.fasta',
        _out = lambda wildcards, output: PurePath(output[0]).parent
    threads: 8
    resources:
        mem_mb = 1500,
        walltime = '4:00'
    shell:
        '''
         RepeatMasker -xsmall -pa $(({threads}/2)) -dir {param._out} -lib {params.library} -qq -no_is {input} 
        '''

rule split_fasta:
    input:
        rep_out = '/cluster/work/pausch/danang/psd/scratch/assembly/{chrom}/{chrom}_rep/{asm}_{chrom}.fa.out',
        fasta = '/cluster/work/pausch/danang/psd/scratch/assembly/{chrom}/{asm}_{chrom}.fa',
        fai = '/cluster/work/pausch/danang/psd/scratch/assembly/{chrom}/{asm}_{chrom}.fa.fai'
    output:
        trimmed = temp('edit_distance2/{chrom}_{asm}.chunks.trimmed.fasta'),
        untrimmed = temp('edit_distance2/{chrom}_{asm}.chunks.untrimmed.fasta')
    params:
        centromere_min_score = 50000,
        telomere_min_score = 1000,
        window = config.get('window',500000),
        chr_size = lambda wildcards, input: int([l.split()[1] for l in open(input.fai)][0])
    resources:
        mem_mb = 500,
        walltime = '30',
        disk_scratch = 1
    run:
        centro_region = subprocess.run(f'awk \'/Satellite\/centr/ {{print $5"\\t"$6"\\t"$7"\\t"$1}}\' {input.rep_out} | bedtools merge -d 250000 -i - -c 4 -o sum | head -n 1',shell=True,capture_output=True).stdout.decode("utf-8")
        telo_region = subprocess.run(f'awk \'/\\(TTAGGG\\)n/||/\\(TAGGGT\\)n/||/\\(AGGGTT\\)n/||/\\(GGGTTA\\)n/||/\\(GGTTAG\\)n/||/\\(GTTAGG\\)n/ {{print $5"\\t"$6"\\t"$7"\\t"$1}}\' {input.rep_out} | bedtools merge -d 1000 -i - -c 4 -o sum | tail -n 1',shell=True,capture_output=True).stdout.decode("utf-8")

        print(f'Regions: {centro_region=} {telo_region=}')

        start = 1 if (not centro_region or int(centro_region.split()[-1]) < params.centromere_min_score) else int(centro_region.split()[2])
        end = params.chr_size if (not telo_region or int(telo_region.split()[-1]) < params.telomere_min_score) else int(telo_region.split()[1])
        
        masked_regions = ''
        if start != 1:
            masked_regions += '\n'.join((f'{wildcards.chrom}_{wildcards.asm}:{i}-{min(start,i+params.window-1)}' for i in range(1,start+1,params.window)))
        if end != params.chr_size:
            if masked_regions:
                masked_regions += '\n'
            masked_regions += '\n'.join((f'{wildcards.chrom}_{wildcards.asm}:{i}-{min(params.chr_size,i+params.window-1)}' for i in range(end,params.chr_size+1,params.window)))

        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            temp.write(masked_regions)
            temp.seek(0)
            subprocess.run(f'samtools faidx -r {temp.name} {input.fasta} > {output.untrimmed}',shell=True)
    
        with tempfile.NamedTemporaryFile(mode='w+t') as temp:
            temp.write('\n'.join((f'{wildcards.chrom}_{wildcards.asm}:{i}-{min(end,i+params.window-1)}' for i in range(start,end+1,params.window))))
            temp.seek(0)
            subprocess.run(f'samtools faidx -r {temp.name} {input.fasta} > {output.trimmed}',shell=True)

def get_threads(wildcards,input):
    if input.size_mb == 0:
        return 1
    if wildcards.trimmed == 'trimmed':
        return 12
    else:
        return 2
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
        return 10000
    else:
        return 40000

#for i in {1..29}; do awk 'BEGIN{OFS=FS="\t";}$1!="L"||$6!="*"{print;}$1=="L"&&$6=="*"{$6="0M";print;}' < cactus_ancestral_${i}.gfa > cactus_${i}.gfa; done
rule graphaligner:
    input:
        gfa = '/cluster/work/pausch/to_share/sv_analysis/final_data/graph/{pangenome}/{pangenome}_{chrom}.gfa',
        fasta = lambda wildcards: 'edit_distance2/{chrom}_{asm}.chunks.{trimmed}.fasta'
    output:
        temp('edit_distance2/{chrom}_{asm}.{pangenome}.{trimmed}.gaf')
    params:
        preset = config.get('X-preset','dbg')
    threads: lambda wildcards,input: get_threads(wildcards,input)
    resources:
        mem_mb = lambda wildcards,input: get_memory(wildcards,input),
        walltime = '2:00'
    shell:
        '''
        if [ -s {input.fasta} ]; then
          GraphAligner -g {input.gfa} -f {input.fasta} --multimap-score-fraction 1 -a {output} -t {threads} --discard-cigar \
          -x dbg -C -1 \
          --max-trace-count 5 --precise-clipping 0.9 \
          --seeds-minimizer-ignore-frequent 0.001 -b 10
        else
          touch {output}
        fi
        '''

rule bin_edit_distance:
    input:
        'edit_distance2/{chrom}_{asm}.{pangenome}.{trimmed}.gaf'
    output:
        'edit_distance2/{chrom}_{asm}.{pangenome}.{trimmed}.dist'
    shell:
        '''
         awk '{{M[$1]+=$10;L[$1]+=$11}} END {{ for (key in M) {{ print key,M[key],L[key] }} }}' {input} | sort -k1,1V > {output}
        '''

#rule defunct
#    shell:
#        '''
#        CENT=$(awk '/Satellite\/centr/ {{print $5"\\t"$6"\\t"$7"\\t"$1}}' {input} | bedtools merge -d 250000 -i - -c 4 -o sum | head -n 1 | cut -f 3,4)
#        TEL=$(awk '/\\(TTAGGG\\)n/||/\\(TAGGGT\\)n/||/\\(AGGGTT\\)n/||/\\(GGGTTA\\)n/||/\\(GGTTAG\\)n/||/\\(GTTAGG\\)n/ {{print $5"\\t"$6"\\t"$7"\\t"$1}}' {input} | bedtools merge -d 1000 -i - -c 4 -o sum | tail -n 1 | cut -f 2,4)
#
#        START=1
#        if [[ ${{CENT#*-}} -gt {params.centromere_min_score ]]
#        then
#          START=${{CENT%-*}}
#        fi
##        END={params.size}
#        if [[ ${{TEL#*-}} -gt {params.telomere_min_score ]]
#        then
#          END=${{TEL%-*}}
#        fi
#        '''
