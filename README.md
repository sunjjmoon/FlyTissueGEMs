# Fly Tissue-specific genome-scale metabolic models (FlyTissueGEMs)

This repository contains the MATLAB and R code used to reconstruct 32 tissue-specific genome-scale metabolic models (GEMs) for *Drosophila melanogaster* and to perform flux balance analysis (FBA), flux variability analysis (FVA), flux sampling, sensitivity analysis, NADH demand analysis, and related visualizations.

## Repository structure
- **1_Gene_EC_annotation_update/**  
  Scripts for mapping gene–enzyme associations and updating EC annotations.

- **2_making_tissue_specific_GEMs/**  
  Scripts for generating tissue-specific GEMs using the updated generic fruitfly model (fruitfly3.mat) and pseudo-bulk transcriptomic data from Fly Cell Atlas (FCA). Includes scripts for subsystem coverage and metabolic task (functional) analysis.

- **3_flux_analysis/**  
  Scripts for FBA, pFBA, FVA, flux sampling, NADH demand, sensitivity analysis, and pathway-level flux analysis.

- **Models/**  
  Base generic model (fruitfly3.mat).

- **Tissue_specific_GEMs/**  
  Reconstructed 32 tissue-specific GEMs (.mat). 

## Software and main packages
- MATLAB R2020b (or newer)
- RAVEN Toolbox (v2.10. or newer, for model reconstruction)
- COBRA Toolbox (v3.0 or newer, for flux analysis)
- R (v4.2.1 or newer)

## Code Availability
The  version of the code used in this study is archived on Zenodo:

**DOI: 10.5281/zenodo.17684286**

## Citation
If you use this code or models, please cite:

*Moon SJ et al., "Modeling tissue-specific Drosophila metabolism identifies high sugar diet-induced metabolic dysregulation in muscle at reaction and pathway levels (2025)" https://doi.org/10.5281/zenodo.17684286*

---
