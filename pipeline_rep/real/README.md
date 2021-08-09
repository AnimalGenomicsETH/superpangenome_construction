
### Whole genome-graph construction with real data

This pipeline will create a whole-genome autosomes graph from multiple assemblies:
The pipelines work as follow:
1. Split genome into a chromosome
2. Create chromosome-level graph, since many graph tools does not yet scale to whole-genome application
    - Additional step for minigraph: it will calculate distance between assemblies and determine order of the graph 
3. Combine into a whole-genome graph
4. Create a graph modification accordingly to make it compatible with vg operations, and create index required for vg mapping and variant calling 

### ! Important 
- *The assembly naming and the contig name specification*

Assembly name: whole genome autosomes `Angus_aut.fa`. Should be capitalized and with underscore. *Tested only on autosome*

Contig name: Should be unique, advised to use the sample identifier in the contig name. For reference backbone, you can use numeric assigned (1-29). 
To add sample identifier you can use `fa_renamer.py`. 


