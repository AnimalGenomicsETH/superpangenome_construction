cat trimmed_regions.csv | tr '\t' '\n' | sed 's/BS1/OBda/g;s/BS2/OBdb/g;s/BS3/OBsa/g;s/BS4/OBsb/g;s/BS5/BS1a/g;s/BS6/BS1b/g;s/BS7/BS2a/g;s/BS8/BS2b/g' | sort -t ',' -o trimmed_regions.csv -k1,1 -k2,2n


#comparing records
for i in {1..29}; do bcftools stats /cluster/work/pausch/to_share/sv_analysis/final_data/snarl/snarl_new/minigraph_gt/minigraph_${i}_snarl_norm_nodup.vcf | awk '$1=="SN"&&/number of records/ {print $6}'; done | awk '{c+=$1} END {print c}'

#alleles per bubble
for i in {1..29}; do bcftools norm -D -Oz /cluster/work/pausch/to_share/sv_analysis/final_data/snarl/snarl_new/minigraph/minigraph_${i}_snarl_norm.vcf  | bcftools query -f '%ID\n' | sed 's/_\w*//g' | sort | uniq -c | awk '{print $1}';done | sort | uniq -c | awk '{print $2,$1,"Default"}'


## Approximately count the number of pathed-alleles in each top-level bubble
## Bubble ID is not always encoded in a consistent manner, so some are top-level and some are nested varation (but top-level within a subbubble)
## Output is in minigraph_alleles.csv

for i in {1..29}; do sed -i -r 's/(LV=[0-9])[,0-9]*;/\1;/' P-line.${i}.normed.vcf; done
echo "Alleles/Bubble Count minigraph" > minigraph_alleles.csv
for p in Default P-line
do
  for i in {1..29}
  do
    sort <(bcftools view -e 'INFO/LV>0' /cluster/work/pausch/alex/GFA_VNTRs/minigraph_vcfs/${p}.${i}.normed.vcf | bcftools query -f '%ID\n' | sed 's/_\S*//g' | sed 's/,\S*//g') \
         <( bcftools view -i 'INFO/LV>0' /cluster/work/pausch/alex/GFA_VNTRs/minigraph_vcfs/${p}.${i}.normed.vcf | bcftools query -f '%INFO/PS\n' | sed 's/_\S*//g' | sed 's/,\S*//g') \
    | uniq -c | awk '{print $1}'
  done | sort | uniq -c | awk -v P=${p} '{print $2,$1,P}' >> minigraph_alleles.csv
done


## Find nodes present in GFA but missing in bed files
## These are paths that are included in the graph, but post-hoc unable to be realigned
## Output in missing_nodes.csv

echo "chromosome sample count" > missing_nodes.csv
for c in {1..29}
do
  awk '$1=="S" {print $2,$5,$4}' /cluster/work/pausch/to_share/sv_analysis/final_data/graph/minigraph/minigraph_${c}.gfa > gafs.txt
  grep -hoP "s\d+" /cluster/work/pausch/to_share/sv_analysis/final_data/snarl/snarl_new/minigraph/call_sv/*/*_${c}.bed | sort -Vu > beds.txt
  GFA=$(wc -l gafs.txt | cut -f 1 -d ' ')
  grep -Ff <(comm -13 <(sort beds.txt) <(cut -f 1 -d' ' gafs.txt | sort)) gafs.txt | grep -v UCD | awk -v c=${c} -v g=$GFA '{split($2,a,"_");split($3,b,":"); print c,a[2],b[3],g}' >> missing_nodes.csv
done


grep -hoP "(SUPP_VEC=\d+|INTRASAMPLE_IDLIST=\S+;)" *_B_merged.vcf | paste -d " "  - - | awk '{split($2,a,"."); printf $1" "; for(i in a) printf gsub(",","X",a[i])+1" "; print ""}'

for j in pggb cactus; do for i in {1..29}; do bcftools norm -m -any ${j}.${i}.vcf | bcftools view -i 'abs(ILEN)>=50' -a | bcftools stats | grep SN >> info.txt ; done; done
awk '/number of records/||/number of no/ {print $6}' info.txt | paste -d " "  - - | awk '{print $2/$1}' | awk 'NR>29{c+=$1;n+=1} END {print c/n*100}'



## Getting the VNTR files


awk 'NR == FNR { a[$0]; next } {if($0 in a) {print $0",yes";next}{print $0",no"}} ' minigraph.trs all.trs > p2.trs

awk 'NR>1 {print $1","$2","$3}' Nellore.naive.bed > ../advntr.trs

grep -Ff ../advntr.trs ../p2.trs > ../p3.trs 

awk 'NR == FNR { a[$0];next } {if($0 in a) {print $0",yes";next}{print $0",no"}} ../p3.trs ../p2.trs > ../p4.trs
