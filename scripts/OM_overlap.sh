bcftools merge -o OM.vcf.gz *gz
bcftools view -o OM.filtered.vcf.gz -i 'ALT="<DEL>"||ALT="<INS>"||ALT="<DUP>"'  OM.vcf.gz
paste <(zgrep -v "#" OM.filtered.vcf.gz | cut -f 1-2 | sed 's/chr//g') <(zgrep -oP ";END=\d*" OM.filtered.vcf.gz | cut -c 6-) | awk '$3>$2 {print $0;next} {print $1"\t"$3"\t"$2}' | sort -k1,1n -k2,2n | bedtools merge -i - > OM.merged.bed
awk '$3-$2<1000000' OM.merged.bed > OM.filtered.bed


bcftools view -i '(ALT="<DEL>"||ALT="<INS>"||ALT="<DUP>")&&ABS(INFO/SVLEN)<1000000'  OM.vcf.gz | bcftools norm -d none -o OM.filtered.vcf.gz
