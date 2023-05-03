[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7891567.svg)](https://doi.org/10.5281/zenodo.7891567) ![Build](https://github.com/AnimalGenomicsETH/superpangenome_construction/actions/workflows/snakemake.yaml/badge.svg)
# Graph construction method impacts variation representation and analyses in a bovine super-pangenome

This work focuses on comparing pangenome constructions methods for the downstream implications of how well (and truthfully) variation is converted into currently-useful formats like VCF.

This reposistory is currently a work-in-progress as code is centralised and cleaned up.

## Components

This project involves multiple steps. We first
- construct pangenomes from 12 bovine assemblies with
  - minigraph
  - pggb
  - cactus
- decompose pangenomes into vcf
- compare vcf files for accuracy to assembly-called _truth_
- assess applications to new topics like VNTRs.


### Citation
The preprint associated with this work can be found [here](https://www.biorxiv.org/content/10.1101/2022.09.17.508368).
