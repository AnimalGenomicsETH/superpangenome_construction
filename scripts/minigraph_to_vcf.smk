rule all:
    input:
        #expand('graph_{c}.L{L}.wave.SV.vcf.gz',c=range(1,30),L=(2,5,10,30,50)),
        expand('jasmine.{chr}.vcf',chr=range(1,30))

rule vg_deconstruct:
    input:
        'graph_{chr}.L{L}.gfa'
    output:
        wave = 'graph_{chr}.L{L}.wave.vcf.gz',
        SV = 'graph_{chr}.L{L}.wave.SV.vcf.gz'
    threads: lambda wildcards: {50:6,30:6,10:8,5:12,2:18}[int(wildcards.L)]
    resources:
        mem_mb = 4000,
        disk_scratch = 25
    shell:
        '''
        vg deconstruct -t {threads} -a -p {wildcards.chr}_UCD {input} |\
        vcfwave -t {threads} |\
        sed 's/\\tFORMAT//g' |\
        bcftools norm --threads {threads} -f /cluster/work/pausch/danang/psd/scratch/assembly/{wildcards.chr}/UCD_{wildcards.chr}.fa -m - |\
        bcftools sort -T $TMPDIR -o {output.wave} -
        tabix -fp vcf {output.wave}
        bcftools stats --threads {threads} {output.wave} > {output.wave}.stats
        bcftools view -i 'abs(ILEN)>=50' -o {output.SV} {output.wave}
        tabix -fp vcf {output.SV}
        bcftools stats --threads {threads} {output.SV} > {output.SV}.stats
        '''

rule jasmine_intersect:
    input:
        expand('graph_{{chr}}.L{L}.wave.SV.vcf',L=(2,5,10,30,50))
    output:
        vcf = 'jasmine.{chr}.vcf',
        #_count = 'jasmine.{chr}.count'
    params:
        _input = lambda wildcards, input: ','.join(input)
    threads: 4
    resources:
        mem_mb = 2000,
        disk_scratch = 10
    conda:
        'jasmine'
    shell:
        '''
        jasmine --comma_filelist file_list={params._input},/cluster/work/pausch/to_share/sv_analysis/snarl/merged_assembly/vcf_split/{wildcards.chr}_assembly_norm.vcf threads={threads} out_file={output.vcf} out_dir=$TMPDIR genome_file=/cluster/work/pausch/danang/psd/scratch/assembly/{wildcards.chr}/UCD_{wildcards.chr}.fa min_seq_id=.5 --pre_normalize --ignore_strand --allow_intrasample --normalize_type
        '''
#grep -vE "SVTYPE=(INV|TRA)" {output.vcf} | grep -oP "(SVLEN=-?\d*|SUPP_VEC=\d{2})" | sed 's/[A-Z,=,_]*//g'  | paste -s -d' \n' > {output._count}
