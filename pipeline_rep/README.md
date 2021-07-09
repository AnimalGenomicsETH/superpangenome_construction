# Compare graphs based on the simulated structural variations

The pipeline will:
- Simulated SV and create a new genome
- Integrate this genome with initial genome into a pggb , minigraph, and cactus graphs
- Deconstruct graph and report variation not found in the reference path
- Compare the variation found with the truth variations from simulation


Tested graphs (Documented 9 July):
- Minigraph (version 0.15-r426)
- PGGB (snapshot docker of `docker pull ghcr.io/pangenome/pggb:20210706131009d4e06f`) 
- Cactus version 2.0.1