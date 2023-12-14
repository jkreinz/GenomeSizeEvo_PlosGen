# GenomeSizeEvo_PlosGen
Scripts and metadata for Genome Size Evolution in Amaranthus Paper:
The scripts in this main repo are for the analysis of TE density data, repeat content, and genome size, as they pertain to differences among populations by geography, environment, and sex.

4 scripts are included here for reproducing key figures and results in the PlosGen2023 paper:
- TEstats_AcrossTheGenome_Fig1.R: Visualizing the distribution of repeats across the genome and collecting statistics of repeat content (by size, distance from nearest gene, recombination rate, coding sequence density, etc)
- RepeatContent_andGenomeSize_byFloweringTime_Analyses_Fig2Fig3.R: Contains many analyses presented with the paper, ranging from repeat abundance and genome size correlations with predictor variables, to flowering time and growth rate correlations with genome size.
- CopyNumberGWAS_Figure4AandB.R: For running a copy number GWA of flowering time, and investigation associations of a large effect CNV with flowering.
- JointGeneticAnalysis_PGS_GS_CNV_Fig4C.R: Calculate polygenic scores from a GWA of flowering time, and compare it to genome size and previously discovered CNV as explanatory variables of flowering time.

Metadata and input files required for these scripts can be found on [Zenodo](10.5281/zenodo.10366486) 
