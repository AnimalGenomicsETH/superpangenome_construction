#!/usr/bin/env bash 

ulimit -u 10000 

input=$1
ref=$2
output=$3
chromo=$4
mode=$5

singularity exec --bind $(dirname $output):/output \
--bind /usr/lib/locale/ \
--bind $TMPDIR:$TMPDIR \
--bind $(dirname $input):/input \
--bind $(dirname $ref):/ref \
/cluster/work/pausch/danang/psd/bin/sif/deepvariant_1.3.0.sif \
/opt/deepvariant/bin/run_deepvariant --model_type $mode \
--ref /ref/$(basename $ref) \
--reads /input/$(basename $input) \
--output_vcf /output/$(basename $output) \
--num_shards 10 \
--regions ${chromo} \
--intermediate_results_dir $TMPDIR \
