
og=$(mktemp).og

echo "odgi extract -i ${1}.gfa -o $og -r HER:${2}-${3} -d 1000 -t 2"

odgi extract -i ${1}.gfa -o $og -r HER:${2}-${3} -d 2500 -t 2

og2="${1}_${2}-${3}.og"
odgi sort -i ${og} -O -t 2 -o ${og2}
for i in $(odgi paths -L -i ${og2})
do
  odgi depth -i ${og2} -d -s <(echo "${i}") > depth.${i}.txt
done

(echo query.name query.start query.end ref.name ref.start ref.end score inv self.cov n.th |\
  tr ' ' '\t'; odgi untangle -i ${og2} -r HER:${2}-${3} --threads 2 -m 100  | bedtools sort -i - ) |\
  awk '$8 == "-" { x=$6; $6=$5; $5=x; } { print }' |\
    awk -v OFS='\t' '{print $1,$4":"$5"-"$6,$2,$3,"1"}' | sed 's/query.*/molecule\tgene\tstart\tend\tstrand/' > ${1}_${2}-${3}.gggenes.txt

odgi view -i ${og2} -t 2 -g > ${1}_${2}-${3}.gfa
