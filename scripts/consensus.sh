#!/usr/bin/env bash

singularity exec --bind $PWD $SIFDIR/pggb.sif \
smoothxg -t 20 -C 1000 -g pggb_100000/25_comb.fa.832b3e4.34ee7b1.seqwish.gfa -C cons,10,100,1000 \
-o pggb_100000/test_smooth.gfa -w 10000000 -K -d 2000 -I 0 -R 0 -j 100 -e 0 -l 5000 -p 1,7,11,2,33,1 -Q Consensus_ -V 

