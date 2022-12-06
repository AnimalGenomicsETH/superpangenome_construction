

echo -e "pangenome\tCPU\tmemory\tsize"
for P in minigraph pggb cactus
do
  MEM=$(grep -l "Successfully" logs/${P}_construct/chromosome-*202212*out | xargs grep -h "Max Memory"| awk '$4>c{c=$4} END {print c/1000}')
  CPU=$(grep -l "Successfully" logs/${P}_construct/chromosome-*202212*out | xargs grep -h "CPU time"| awk '{c+=$4} END {print c/3600}')

  SIZE=$(du -c graphs/${P}/{1..29}.gfa | awk '/total/ {print $1/1048576}')
  echo -e "$P\t$CPU\t$MEM\t$SIZE"
done
grep -l "Successfully" logs/repeatmasker_soft/chromosome-*202212*out | xargs grep -h "CPU time"| awk '{c+=$4} END {print c/3600}'
grep -l "Successfully" logs/cactus_convert/chromosome-*202212*out | xargs grep -h "CPU time"| awk '{c+=$4} END {print c/3600}'
