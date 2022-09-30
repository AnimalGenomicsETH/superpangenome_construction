## LOCAL ASSEMBLIES
for ID in GAU,GxP.hap1 PIE,GxP.hap2 NEL,NxB.hap1 BSW,NxB.hap2 OBV,OxO.hap2
do
  seqtk subseq /nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/assemblies/CCS/${ID#*,}.hifiasm.fasta <(echo {1..29} | tr ' ' '\n') > ${ID%,*}.fasta
done

##LOCAL COPY OF REFERENCE
seqtk subseq /cluster/work/pausch/inputs/ref/BTA/UCD1.2/ARS-UCD1.2_Btau5.0.1Y.fa <(echo {1..29} | tr ' ' '\n') > ref.fasta


## OTHER ASSEMBLIES AVAIBLE VIA NCBI
for ID in BRA,GCF_003369695.1 ANG,GCA_003369685.2 SIM,GCA_018282465.1 HIG,GCA_009493655.1 YAK,GCA_009493645.1 BIS,GCA_018282365.1
do
  curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v1/genome/accession/${ID#*,}/download?filename=${ID#*,}.zip" -H "Accept: application/zip"
  unzip -o ${ID#*,}.zip && rm ${ID#*,}.zip
  for i in {1..29}; do cat ncbi_dataset/data/${ID#*,}/chr${i}.fna | sed "s/>.*/>${i}/g" | seqtk seq >> ${ID%,*}.fasta; done
done

rm -rf ncbi_dataset



#for p in BRA ANG
#do
#  sort -k11,11nr ${p}.paf | grep -vE "(NKL|X|Y)" | sort -k6,6n | awk '{print $6,$5,$11}' | sort -k2,2h | awk '{seen[$1" "$2]+=$3} END { for (key in seen) { print key,seen[key] } }' |  sort -k1,1n -k3,3nr | awk '!orient[$1]{orient[$1]=$2} END { for (key in orient) { print key,orient[key] } }' > ${p}.directions
#  awk '$2=="-" {print $1}' ${p}.directions | seqtk subseq ${p}.fasta - | seqtk seq -r - > ${p}_fixed.fasta
#  awk '$2=="+" {print $1}' ${p}.directions | seqtk subseq ${p}.fasta - >> ${p}_fixed.fasta
#rm ${p}.directions
#mv ${p}_fixed.fasta ${p}.fasta
#done
