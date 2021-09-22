# Comparison of graph-building approaches based on the simulated structural variations

The genome graph built from three tools (pggb, cactus, and minigraph) are compared in three ways
1. `run_replicate.smk`: Deconstruct variations directly from graphs (paths in bubbles not part of reference nodes) 
2. `map_call.smk`: Map Illumina-short read data and call variations from the alignments
3. `map_lr.smk`: Map long-read (Pacbio and ultra long ONT) and call variations from alignments


First the pipeline in `run_replicate.smk` will:
- Simulate SV and create a new genome
- Integrate this genome with initial genome into a pggb , minigraph, and cactus graphs
- Deconstruct graph and report variation not found in the reference path
- Compare the variation found with the truth variations from simulation

In `map_call.smk` short-read data simulated with `Mason` are mapped to cactus, pggb, and minigraph graph and SV called with `vg call`. For comparison short-read were also mapped to linear genome with `bwa` and SV called with `delly` and `manta`. To differentiate the improvement from mapping difference the pipeline will also map to primary graph containing only reference backbone. 

In `map_lr.smk` long-read data simulated with `pbsim2` of both Pacbio and ONT using model based on P6C4.model and R95 chemistry respectively. Simulated reads mapped to three graphs with `GraphAligner` and variations were called with `vg call`. To compare with linear mapping, the pipeline mapped simulated reads to the reference backbone with  `ngmlr` and call SV with `sniffles`. 

### Tested graphs (as 9 July):
- Minigraph (version 0.15-r426)
- PGGB (snapshot docker of `docker pull ghcr.io/pangenome/pggb:20210706131009d4e06f`) 
- Cactus version 2.0.1

### Tested tools (as 21 July)
- vg: v1.33.0 "Moscona"
- Mason, mason_simulator module  2.0.9 [e165baf]
- BWA 0.7.17-r1188
- delly  0.8.7
- manta 1.6.0
- pbsim2 e71f789
- GraphAligner 1.0.13
- ngmlr 0.2.7 
- sniffles 1.0.12
- SURVIVOR 1.0.7

