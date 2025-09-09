# Scripts

This folder contains R scripts to reproduce all analyses and figures for the manuscript.  
Scripts are numbered to match the figure panels in the paper. Each script reads from `../data/` and writes to `../figures/`.

## Contents

- `01_tailocin_killing_matrix.R`  
  Plot phylogenetic tree + killing matrix heatmap.

- `02_HTF_OBC_sensitivity_summary.R`  
  Summarize killing assay results by OBC presence/absence; Wilcoxon tests.

- `03_plant_health_OBC.R`  
  Compare mean green pixels (plant health) between OBC groups.

- `04_tnseq_scatter.R`  
  Visualize TnSeq-derived plant fitness as a scatterplot.

- `05_geneplot.R`  
  Plot gene map of the focal OBC cluster colored by TnSeq results.

- `06_plant_infections_luciferase.R`  
  Analyze bacterial load (luminescence fold-change) at 7 dpi; mixed models and pairwise contrasts vs WT.

- `07_plant_disease.R`  
  Analyze qualitative plant outcomes (healthy vs diseased/dead); logistic regressions and Dunnett-style contrasts vs WT.

- `08_disease_vs_load.R`  
  Relate plant health proportions to bacterial load; per-ecotype linear models and correlation plot.

- `09_oxidative_stress.R`  
  Placeholder script for oxidative stress survival assays (in progress).

