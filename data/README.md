# Data directory

This directory contains all raw input data and processed output files used in the analyses
for the Ps1524 tailocin trade-offs manuscript.

All data are organized into three top-level folders:

- `input/` – raw or minimally processed data used as inputs to analysis scripts  
- `output/` – processed data products, summary tables, and statistical results  
- `p25.C2_references/` – reference genome files for the focal strain p25.C2

File usage follows the manuscript and script order.

---

## input/

### Tailocin killing assays (Fig. 1A–B)
- `tailocin_killing_matrix.csv` – pairwise tailocin killing assay results
- `tailocin_killing_long.csv` – long-format killing data for plotting/statistics
- `HTF_length_haplotype.csv` – tail fiber length haplotype assignments
- `ps_1524_uncollapsed_5_2018.nwk` – phylogenetic tree used in killing matrix plots

---

### Plant health and infection assays (Fig. 1C, Fig. 3C, Fig. S2–S4)
- `plant_health_OBC_mean_pixels.csv` – mean green pixel counts per strain
- `plant_health_OBC_presence.csv` – O-antigen biosynthesis cluster presence/absence
- `plant_infections_luciferase_7dpi.csv` – luciferase fold-change at 7 dpi
- `plant_infections_layout.csv` – plate layout metadata for luciferase assays
- `plant_infections_contam.csv` – contamination flags
- `plant_disease.csv` – qualitative disease scoring
- `plant_CFUs_Synthesis.csv` – synthesized CFU counts across ecotypes, doses, and timepoints

---

### TnSeq analyses (Fig. 2B–C)
- `tnseq_plant.csv` – in planta TnSeq read counts
- `tnseq_tailocin.csv` – tailocin TnSeq fitness data
- `tnseq_ortho.csv` – ortholog mappings between p25.C2 and DC3000
- `geneplot_p25C2.csv` – O-antigen biosynthesis cluster gene coordinates and annotations

---

### In vitro growth curves (Fig. S4)
- `invitro_growth_synthesis.csv` – OD600 time-series data for WT and mutants

---

### Aggregation phenotypes (Fig. S5)
- `aggregation.csv` – aggregation phenotype calls (R1–R3)

---

### Antimicrobial peptide assays (Fig. S6)
- `AMP_synthesis.csv` – disk diffusion inhibition area measurements across strains, AMPs, and concentrations

---

### Oxidative stress / ROS assays (Fig. S7)
- `ROS_layout_R12.csv` – plate layout metadata
- `ROS_layout_R13.csv`
- `ROS_layout_R14.csv`
- `ROS_layout_R15.csv`
- `ROS_layout_R16.csv`
- `ROS_OD600_R12.csv` – OD600 growth curves
- `ROS_OD600_R13.csv`
- `ROS_OD600_R14.csv`
- `ROS_OD600_R15.csv`
- `ROS_OD600_R16.csv`

---

### qPCR raw data (Fig. S8)
- `qPCR_RFUs.csv` – raw fluorescence (RFU) values
- `qPCR_Template.csv` – plate templates and sample metadata
- `qPCR_SynthesisMNE.csv` – merged, plate-level processed qPCR data prior to final summarization

---

## output/

### Tailocin killing summary (Fig. 1B)
- `tailocin_killing_summary_stats.txt` – statistical comparisons of sensitivity by OBC status

---

### Plant health (Fig. 1C)
- `plant_health_OBC_wilcoxon_summary.txt` – Wilcoxon test results comparing plant health across OBC groups

---

### TnSeq trade-off analyses (Fig. 2B)
- `tnseq_merged.csv` – merged plant × tailocin gene-level fitness dataset
- `tnseq_manuscript_table.csv` – list of significant resistance trade-off genes

---

### In vitro growth curves (Fig. S4)
- `invitro_growthcurveR_data.csv` – Growthcurver-derived parameters (r, K, AUC)
- `invitro_growthcurve_stats.csv` – statistical comparisons vs WT

---

### Aggregation phenotypes (Fig. S5)
- `aggregation_stats.csv` – pairwise Fisher’s exact tests vs WT
- `aggregation_chisq_overall.csv` – overall χ² test across strains

---

### Antimicrobial peptide assays (Fig. S6)
- `AMP_analysis.csv` – summarized inhibition areas and statistics
- `AMP_PMB_stats.csv` – polymyxin B–specific statistical results
- `AMP_PME_stats.csv` – polymyxin E–specific statistical results

---

### Oxidative stress / ROS assays (Fig. S7)
- `ROS_WT_MIC_by_replicate.csv` – per-replicate WT MIC estimates
- `ROS_MIC_table.csv` – strain-level MIC summaries
- `ROS_emmeans_contrasts_vsWT_scaledMIC.csv` – mixed-model contrasts vs WT

---

### qPCR processed outputs (Fig. S8)
- `qPCR_Cq_Eff.csv` – raw Cq and efficiency estimates
- `qPCR_Cq_Eff_cleanup.csv` – filtered Cq/efficiency table after QC
- `qPCR_checks.csv` – quality control diagnostics
- `qPCR_NE_Mock.csv` – normalized expression values (mock)
- `qPCR_NE_Time0.csv` – normalized expression values (time 0)
- `qPCR_MNE_summary_Mock.csv` – summarized MNE values (mock)
- `qPCR_MNE_summary_Time0.csv` – summarized MNE values (time 0)
- `qPCR_SynthesisMNE_summary.csv` – final summarized qPCR dataset used for plotting and statistics

---

## p25.C2_references/

Reference genome files for the focal strain *Pseudomonas syringae* p25.C2.

- `p25.C2.fna`
- `p25.C2.faa`
- `p25.C2.gff`
- `p25.C2.gbk`
- `plate25.C2.pilon.contigs_renamed.fasta`

