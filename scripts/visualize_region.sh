#ususr/bin/env bash

usage (){

cat << EOF
##################################################
        visualize the 
        -p prefix
	-c chr:pos1[-pos2]
	-e expand context by x node
        # Note that this required vg v1.34, Bandage, and montage for combining
        # can be activated with `conda activate imageman` 
##################################################
EOF

exit 1

}

if [[ $# -eq 0 ]];then usage; fi

while getopts ":p:c:e:" opt;do
        case $opt in
                p) prefix=$OPTARG ;;
		c) chromo=$OPTARG ;;
		e) context=$OPTARG ;;
                \?|:|h) usage;;
        esac
done

shift $((OPTIND - 1))


context=${context:-0}
echo $context
contig=$(cut -f1 -d":" <<< $chromo)
pos=$(cut -f2 -d":" <<< $chromo)
 
for grtype in pggb cactus minigraph 
do 
  infile="graph/${grtype}/graph_${grtype}_combined.gfa"
  # subset the graph
  if [[ $grtype == "cactus" ]];then 
	  vg find -p ${contig}.${contig}:${pos} -c ${context} -x ${infile} > ${prefix}_${grtype}.gfa
  else

	  vg find -p ${chromo} -c ${context} -x ${infile} > ${prefix}_${grtype}.gfa
  fi
  # make fasta of the gfa
  awk '$1 ~ /S/{print ">"$2"\n"$3}' ${prefix}_${grtype}.gfa > ${prefix}_${grtype}.fa
  # visualize it 
  Bandage image ${prefix}_${grtype}.gfa ${prefix}_${grtype}.png --names --fontsize 6
  # node label information
  ./node_labels.py -i ${prefix}_${grtype}.gfa -o ${prefix}_${grtype}_label.tsv
done


# Combining images
montage -geometry 900x -tile 3x1 \
   ${prefix}_pggb.png \
   ${prefix}_cactus.png \
   ${prefix}_minigraph.png \
   ${prefix}_comb.png 
