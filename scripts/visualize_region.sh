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

 
for grtype in pggb cactus minigraph 
do 
  infile="graph/${grtype}/graph_${grtype}_combined.gfa"
  # subset the graph
  vg find -p ${chromo} -c ${context} -x ${infile} > ${prefix}_${grtype}.gfa
  # make fasta of the gfa
  awk '$1 ~ /S/{print ">"$2"\n"$3}' ${prefix}_${grtype}.gfa > ${prefix}_${grtype}.fa
  # visualize it 
  Bandage image ${prefix}_${grtype}.gfa ${prefix}_${grtype}.png
done


# Combining images
montage -geometry 900x -tile 3x1 \
   ${prefix}_pggb.png \
   ${prefix}_cactus.png \
   ${prefix}_minigraph.png \
   ${prefix}_comb.png 
