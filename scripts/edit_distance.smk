

import sys
import subprocess
centromere_min_score = 50000
telomere_min_score = 1000
        
input = sys.argv[1]
#'/cluster/work/pausch/alex/assembly/BSWCHEM120151536851/hifiasmv152_100/hap2_split_chrm/5.chrm.fa.out'
centro_region = subprocess.run(f'awk \'/Satellite\/centr/ {{print $5"\\t"$6"\\t"$7"\\t"$1}}\' {input} | bedtools merge -d 250000 -i - -c 4 -o sum | head -n 1',shell=True,capture_output=True).stdout.decode("utf-8")
telo_region = subprocess.run(f'awk \'/\\(TTAGGG\\)n/||/\\(TAGGGT\\)n/||/\\(AGGGTT\\)n/||/\\(GGGTTA\\)n/||/\\(GGTTAG\\)n/||/\\(GTTAGG\\)n/ {{print $5"\\t"$6"\\t"$7"\\t"$1}}\' {input} | bedtools merge -d 1000 -i - -c 4 -o sum | tail -n 1',shell=True,capture_output=True).stdout.decode("utf-8")

        
start = 1 if int(centro_region.split()[-1]) < centromere_min_score else centro_region.split()[2]
end = 9999999999 if int(telo_region.split()[-1]) < telomere_min_score else telo_region.split()[1]
print(f'asm:{start}-{end}')

