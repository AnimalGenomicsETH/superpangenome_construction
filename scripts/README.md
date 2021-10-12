
## Notes on the script usage 


To visualize specific regions on the graph use as follow (an example):

```

bsub -J viz -XF -n 10 -W 01:00 -o viz.log -R "rusage[mem=2000]" \
./visualize_region.sh -p vizreg/test_kit -c 6_UCD:70000000 -e 6

```

p is the prefix of output files    
c is the chromosomal region to visualize     
e indicates how many nodes extended from the regions for visualization    

*Note*  The batch submission should include `-R XF` for X11 forwarding, required by Bandage
It will output:     
-subset graph (in gfa)       
-visualize graph from three program side by side (in comb file)          
-label of each nodes in `*labels.tsv`    


