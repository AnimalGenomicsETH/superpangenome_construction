#!/usr/bin/env bash

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
chrom_id=$(cut -f1 -d"_" <<< $contig)
pos=$(cut -f2 -d":" <<< $chromo)
 
for grtype in pggb cactus minigraph 
do 
  infile="graph/${grtype}/graph_${grtype}_combined.gfa"
  # subset the graph
  if [[ $grtype == "cactus" ]];then 
	  infile="graph/cactus/${chrom_id}/cactus_drop_${chrom_id}.gfa"
	  vg find -p ${contig}.${contig}:${pos} -c ${context} -x ${infile} > ${prefix}_${grtype}.gfa
  elif [[ $grtype == "pggb" ]];then 
	  infile="graph/pggb_${chrom_id}/pggb_${chrom_id}.gfa"
	  vg find -p ${chromo} -c ${context} -x ${infile} > ${prefix}_${grtype}.gfa
  elif [[ $grtype == "minigraph" ]];then 
	  infile="graph/minigraph/graph_${chrom_id}.gfa"
	  vg find -p ${chromo} -c ${context} -x ${infile} > ${prefix}_${grtype}.gfa
  fi
  # make fasta of the gfa
  awk '$1 ~ /S/{print ">"$2"\n"$3}' ${prefix}_${grtype}.gfa > ${prefix}_${grtype}.fa
  # visualize it 
  Bandage image ${prefix}_${grtype}.gfa ${prefix}_${grtype}.png --names --fontsize 6
  # node label information
   if [[ $grtype == "minigraph" ]];then

          awk 'NR==FNR && $1 ~ />/{var=substr($0,2,length($0));arr["s"var];next}$1 in arr{print $0}' \
          ${prefix}_${grtype}.fa graph/minigraph/${chrom_id}/${chrom_id}_comb_coverage.tsv > ${prefix}_${grtype}.tmp

          cat <( echo node_id node_len node_label short_label ref_status ) \
          <( awk '{print $1,$2,$(NF-2),$(NF-1),$NF}' OFS="\t" ${prefix}_${grtype}.tmp) > ${prefix}_${grtype}_label.tsv

          if [[ -f ${prefix}_${grtype}.tmp ]];then rm ${prefix}_${grtype}.tmp; fi

  else
        ./node_labels.py -i ${prefix}_${grtype}.gfa -o ${prefix}_${grtype}_label.tsv
  fi
  #./node_labels.py -i ${prefix}_${grtype}.gfa -o ${prefix}_${grtype}_label.tsv
done


# Combining images
montage -geometry 900x -tile 3x1 \
   ${prefix}_pggb.png \
   ${prefix}_cactus.png \
   ${prefix}_minigraph.png \
   ${prefix}_comb.png 
