#!/usr/bin/env bash 


chromo=$1 
prefix=/cluster/work/pausch/danang/psd/scratch/assembly/${chromo}

vg deconstruct -p UCD.${chromo}_UCD -t 20 -a -e cactus/cactus_simple_${chromo}.gfa > valid_snps/cactus${chromo}_var_graph.vcf 

awk '$11 !~ "0" ' valid_snps/cactus${chromo}_var_graph.vcf > bsw_only_cactus_${chromo}.tsv 

vg deconstruct -p ${chromo}_UCD -t 20 -a -e pggb_${chromo}/pggb_${chromo}.gfa > valid_snps/pggb${chromo}_var_graph.vcf  

awk '$11 !~ "0" ' valid_snps/pggb${chromo}_var_graph.vcf > bsw_only_pggb_${chromo}.tsv 

# variations from assembly vs assembly comparison
minimap2 -cx asm5 -t20 --cs ${prefix}/UCD_${chromo}.fa ${prefix}/BSW_${chromo}.fa > valid_snps/UCD_BSW_${chromo}.paf

sort -k6,6 -k8,8n valid_snps/UCD_BSW_${chromo}.paf > valid_snps/UCD_BSW_${chromo}_sort.paf && rm valid_snps/UCD_BSW_${chromo}.paf

paftools.js call valid_snps/UCD_BSW_${chromo}_sort.paf > valid_snps/BSW_${chromo}_varassemb.tsv  

awk 'NR==FNR{arr[$3]=$0;arr[$4]=$0;next}$2 in arr{print $0,arr[$2]}' valid_snps/BSW_${chromo}_varassemb.tsv valid_snps/pggb${chromo}_var_graph.vcf > valid_snps/pggb${chromo}_overlap.tsv

awk 'NR==FNR{arr[$3]=$0;arr[$4]=$0;next}$2 in arr{print $0,arr[$2]}' valid_snps/BSW_${chromo}_varassemb.tsv valid_snps/cactus${chromo}_var_graph.vcf > valid_snps/cactus${chromo}_overlap.tsv
