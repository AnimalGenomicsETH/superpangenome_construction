echo "pangenome chromosome deletions insertions other"
for p in pggb cactus
do
  for c in {1..29}
  do
    echo -n "${p} ${c} "
    echo $(bcftools stats /cluster/work/pausch/to_share/sv_analysis/final_data/snarl/snarl_new/${p}/${p}_${c}_snarl_norm.vcf |\
    awk '{if($5~/others/){O+=$6}else{if($1=="IDD"){if($3~/-/){D+=$4}else{I+=$4}}}} END {print D,I,O}')
  done
done

for c in {1..29}
do
  echo -n "minigraph ${c} "
  echo $(bcftools stats /cluster/work/pausch/to_share/sv_analysis/final_data/snarl/snarl_new/minigraph_gt/minigraph_${c}_snarlgt_norm.vcf |\
  awk '{if($5~/others/){O+=$6}else{if($1=="IDD"){if($3~/-/){D+=$4}else{I+=$4}}}} END {print D,I,O}')
done
