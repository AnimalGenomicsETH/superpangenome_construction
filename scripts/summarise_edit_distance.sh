echo "pangenome breed chromosome region matched aligned seeded size memory CPU errors"
for p in minigraph cactus pggb
do
  for c in {1..29}
   do
     for b in Angus Bison Brahman BSW Gaur Highland Nellore OBV Pied Simmental UCD Yak BS1 BS2 BS3 BS4 BS5 BS6 BS7 BS8
     do
       for T in trimmed untrimmed
       do 
         echo -n "$p $b $c $T "
         if [ -s edit_distance/${c}_${b}.${p}.${T}.dist ]; then
           echo -n  $(awk '{{split($1,a,":");split(a[2],b,"-"); D+=(b[2]-b[1]);L+=$2;M+=$3}} END {{print L,M,D}}' edit_distance/${c}_${b}.${p}.${T}.dist)
         else
           echo -n "nan nan nan"
         fi
        
         echo " "$(awk -v c=0 '!/>/ {{c+=length($1)}} END {{print c}}' edit_distance/${c}_${b}.chunks.${T}.fasta) \
         $(awk '/Max Memory/ {{print $4}}' logs/graphaligner/chrom-${c}.asm-${b}.pangenome-${p}.trimmed-${T}_*out) \
         $(awk '/CPU/ {{print $4}}' logs/graphaligner/chrom-${c}.asm-${b}.pangenome-${p}.trimmed-${T}_*out) \
         $(grep -hcF Seed logs/graphaligner/chrom-${c}.asm-${b}.pangenome-${p}.trimmed-${T}_*err)
      done
    done
  done
done
