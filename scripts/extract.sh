#!/usr/bin/env bash 

echo "extract path range input out_prefix"

odgi extract -r $1 -L $2 -t 10 -i $3 -o ${4}.og 

odgi view -i ${4}.og -g > ${4}.gfa && rm ${4}.og 

node_labels.py -i ${4}.gfa -o ${4}_label.csv

