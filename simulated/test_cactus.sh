#!/usr/bin/env bash


source /cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/venv/bin/activate
export PATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin:$PATH
export PYTHONPATH=/cluster/work/pausch/danang/psd/bin/cactus-bin-v2.0.1/bin/lib:$PYTHONPATH
chromo=$1


echo -e "${chromo}_OBV\\t_OBV.fa\\n25_sim\\t{input.fasta}" > {params.prefix}/seqfile_{wildcards.rep}.tsv 

cactus-graphmap {params.scrdir}  \
{params.prefix}/seqfile_{wildcards.rep}.tsv \
{input.graph} {output} \
--outputFasta {params.prefix}/cac_{wildcards.rep}.fa \
--realTimeLogging