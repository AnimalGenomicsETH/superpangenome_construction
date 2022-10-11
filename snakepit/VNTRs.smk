

rule odgi_position:
    input:
        gfa = 'graphs/{pangenome}/{chromosome}.gfa',
        TRs = 'VNTRs/TRF/{chromosome}.TRs.bed'
    output:
        'VNTRs/{pangenome}/{chromosome}.liftover.bed'
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '4:00'
    shell:
        '''        
        odgi position -i {input.gfa} -b {input.TRs} -t {threads} > {output}
        '''

