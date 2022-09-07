cat trimmed_regions.csv | tr '\t' '\n' | sed 's/BS1/OBda/g;s/BS2/OBdb/g;s/BS3/OBsa/g;s/BS4/OBsb/g;s/BS5/BS1a/g;s/BS6/BS1b/g;s/BS7/BS2a/g;s/BS8/BS2b/g' | sort -t ',' -o trimmed_regions.csv -k1,1 -k2,2n


#comparing records
for i in {1..29}; do bcftools stats /cluster/work/pausch/to_share/sv_analysis/final_data/snarl/snarl_new/minigraph_gt/minigraph_${i}_snarl_norm_nodup.vcf | awk '$1=="SN"&&/number of records/ {print $6}'; done | awk '{c+=$1} END {print c}'

#alleles per bubble
for i in {1..29}; do bcftools norm -D -Oz /cluster/work/pausch/to_share/sv_analysis/final_data/snarl/snarl_new/minigraph/minigraph_${i}_snarl_norm.vcf  | bcftools query -f '%ID\n' | sed 's/_\w*//g' | sort | uniq -c | awk '{print $1}';done | sort | uniq -c | awk '{print $2,$1,"Default"}'

for i in {1..29}; do sed -i -r 's/(LV=[0-9])[,0-9]*;/\1;/' P-line.${i}.normed.vcf; done
echo "Alleles/Bubble Count minigraph"
for p in Default P-line
do
  for i in {1..29}
  do
    sort <(bcftools view -e 'INFO/LV>0' /cluster/work/pausch/alex/GFA_VNTRs/minigraph_vcfs/${p}.${i}.normed.vcf | bcftools query -f '%ID\n' | sed 's/_\S*//g' | sed 's/,\S*//g') \
         <( bcftools view -i 'INFO/LV>0' /cluster/work/pausch/alex/GFA_VNTRs/minigraph_vcfs/${p}.${i}.normed.vcf | bcftools query -f '%INFO/PS\n' | sed 's/_\S*//g' | sed 's/,\S*//g') \
    | uniq -c | awk '{print $1}'
  done | sort | uniq -c | awk -v P=${p} '{print $2,$1,P}'
done
