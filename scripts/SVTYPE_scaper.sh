echo "pangenome chromosome deletions insertions other"
for p in pggb cactus minigraph assembly
do
  for c in {1..29}
  do
    echo -n "${p} ${c} "
    echo $(bcftools stats ${p}/${c}.vcf |\
    awk '{if($5~/others/){O+=$6}else{if($1=="IDD"){if($3~/-/){D+=$4}else{I+=$4}}}} END {print D,I,O}')
  done
done
