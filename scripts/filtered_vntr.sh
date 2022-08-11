awk -F, -v OFS=, 'NR>1{$4 = sprintf("%.0f", $4); print $1"\t"$4"\t"$4+100000}' /cluster/work/pausch/to_share/sv_analysis/snarl/merged/filtered_excess_variants.csv | sort -k1,1n -k2,2n -u | bedtools slop -i - -g <(cut -f -2 /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai) -b 0 > retained_regions.bed
bedtools complement -i retained_regions.bed -g <(cut -f -2 /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai) | grep -P "^\d" > filtered_regions.bed


bedtools intersect -a filtered_regions.bed -b /cluster/work/pausch/alex/GFA_VNTRs/TR/d0_m100.L100.bed | wc -l > vntr_hits.txt

bedtools makewindows -g <(cut -f -2 /cluster/work/pausch/alex/REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai | grep -P "^\d") -w 100000 | sort -k1,1n -k2,2n > windows.bed
for i  in {1..10000}
do
  shuf -n 8678 windows.bed | sort -k1,1n -k2,2n | bedtools intersect -a - -b /cluster/work/pausch/alex/GFA_VNTRs/TR/d0_m100.L100.bed | wc -l >> vntr_hits.txt
done



#import scipy.stats as ss
#import numpy as np
#import matplotlib.pyplot as plt
#x=np.linspace(6000,10000,4001)
#real,*simulations=[int(i) for i in open('Danang_graphs/vntr_hits.txt')]
#plt.hist(simulations,density=1,bins=100)
#p=ss.norm.fit(simulations)
#plt.plot(x,ss.norm(*p).pdf(x),'r--')
