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
  
- `04_tnseq_plant.R`  
  Process raw in planta TnSeq read counts from *Pseudomonas syringae* p25C2 grown in Col-0 and Eyach 1.5-2 ecotypes.  
  Performs DESeq2 differential abundance analysis, filters low-read samples, maps DC3000 orthologs to p25C2,  
  and merges host fitness results with tailocin TnSeq data to produce a unified dataset for downstream trade-off analyses.  
  Output: `merged_tnseq.csv` in `data/04_tnseq_scatter/output/`.

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
  Quantify tolerance of WT and OPS mutants to oxidative stress (H₂O₂).  
  Calculates per-well AUCs, normalizes by 0 mM controls, scales doses by WT MIC,  
  and fits mixed-effects models with emmeans contrasts vs WT.  
  Outputs: `WT_MIC_by_replicate.csv`, `emmeans_contrasts_vsWT_scaledMIC.csv`, and `09_oxidative_stress.png`.

- `10_aggregation.R`  
  Quantify strain-level aggregation phenotypes (R1–R3).  
  Generates stacked bar plots of phenotype proportions per strain, performs χ² test across strains,  
  and Fisher’s exact tests vs WT with FDR correction.  
  Outputs: `10_aggregation_stats.csv` and `10_aggregation.png`.

- `11_qPCR.R` 
  - `11a_qPCR_processing_example.R`: Example script demonstrating per-plate qPCR processing from raw fluorescence (RFU) data.
    Includes baseline correction, per-well efficiency estimation, quality control of technical replicates, exclusion of failed wells, reference-gene normalization, and calculation of MNE fold-change values.
    Serves as a template for processing all individual qPCR plates prior to synthesis.
  - `11b_qPCR.R`
    Aggregate analysis of synthesized qPCR expression data across experiments.
    Summarizes MNE-normalized fold-change values across biological replicates, generates time-course and treatment-specific visualizations, and performs adaptive statistical testing (ANOVA/Tukey, Welch/Games–Howell, or Kruskal–Wallis/Dunn) for each gene and timepoint.
    Outputs: `SynthesisMNE_summary.csv` and `11_qPCR.pdf`.
  
- `12_invitro_growth.R`
  Analyze in vitro growth curves of wild-type and O-antigen mutants grown overnight at an initial OD₆₀₀ = 0.01.
  Fits growth parameters using Growthcurver, exports rate (r), carrying capacity (k), and AUC values, and performs statistical comparisons (Tukey, Games–Howell, or Dunn tests) between mutants and WT.
  Outputs: `growthcurveR_data.csv`, `growthcurve_stats.csv`, and `12_growthcurves.pdf`.

- `13_plant_CFUs.R`
  Analyze in planta bacterial growth measured as CFU per leaf following infection with WT and O-antigen mutants.
  Generates ecotype- and timepoint-specific boxplots, applies adaptive statistical testing (t-test, Welch t-test, or Wilcoxon rank-sum) to compare WT vs mutant across inoculum concentrations, and exports a unified results table.
  Outputs: `13_plant_CFUs.pdf` and `13_CFU_stats.csv`.