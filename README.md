# Generation of consensus metabolic reconstructions and community-dependent gap-filling using COMMIT for communites sampled from _Arabidopsis thaliana_.

[![DOI](https://zenodo.org/badge/363932874.svg)](https://zenodo.org/badge/latestdoi/363932874)

## Description
This repository contains all files used for community-dependent gap-filling and subsequent analysis
of consensus genome-scale metabolic reconstructions for two soil communities (Bulgarelli et al. 2012, Schlaeppi et al. 2014).

The consensus reconstructions were generated by merging draft reconstructions from four different approaches:
- [KBase](https://www.kbase.us/) (Arkin et al. 2018)
- [RAVEN 2.0](https://github.com/SysBioChalmers/RAVEN) (Wang et al. 2018)
- [CarveMe](https://github.com/cdanielmachado/carveme) (Machado et al. 2018)
- [AuReMe](http://aureme.genouest.org/) (Aite et al. 2018) / [Pathway Tools](http://pathwaytools.com/) (Karp et al. 2016)

## Requirements for COMMIT
- Matlab (tested with versions R2017b and R2020b)
- [CPLEX solver](https://www.ibm.com/analytics/cplex-optimizer) (tested with v12.9)
- R packages: pheatmap, wesanderson, scales, plotrix, igraph, ape
- [COBRA toolbox v3.0](https://github.com/opencobra/cobratoolbox)
- [HMMER](http://hmmer.org/download.html) (tested with v3.2.1)

## Usage
- all settings can be changed in _COMMIT/code/matlab/options.m_
- the script _COMMIT/code/matlab/commit.m_ contains the workflow for the publication
- the script _COMMIT/code/matlab/gap_filling/run_iterative_gap_filling.m_ performs the conditional gap filling procedure of COMMIT
- the script _COMMIT/code/matlab/gap_filling/gap_fill_individual_models.m_ performs gap-filling on the individual reconstructions without considering the community composition

## Reference
Wendering P, Nikoloski Z (2022) COMMIT: Consideration of metabolite leakage and community composition improves microbial community reconstructions. PLOS Computational Biology 18(3): e1009906. https://doi.org/10.1371/journal.pcbi.1009906
