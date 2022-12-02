import sys
import subprocess

for line in sys.stdin:
    *TR_information, bubbles = line.split()
    source, *_, sink = bubbles.split(',')
    print('\t'.join(TR_information),end='\t')

    TR_coordinates = dict()
    
    bedfiles = ' '.join(sys.argv[1:])
    for node,key in zip((source,sink),(4,5)):
        for hit in subprocess.run(f"awk '${key}==\">{node}\"' {bedfiles}",shell=True,capture_output=True).stdout.decode("utf-8").split('\n')[:-1]:
            hit = hit.split()[-1]
            if hit == '.':
                continue
            parts = hit.split(':')
            if parts[3] in TR_coordinates:
                TR_coordinates[parts[3]] += f'-{parts[key]},{parts[2]}'
            elif key == 4:
                TR_coordinates[parts[3]] = parts[key]
            
    print('\t'.join((f'{sample}:{coord}' for sample,coord in TR_coordinates.items())))
