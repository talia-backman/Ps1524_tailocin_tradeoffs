# qPCR plate-level processing (example)
#
# Purpose:
#   This script documents the raw processing pipeline used for all qPCR plates
#   in this study, from fluorescence amplification curves to efficiency-corrected
#   MNE fold-change values.
#
# Scope:
#   This file is a representative example. Identical logic was applied to all
#   qPCR plates, with only plate-specific metadata (targets, wells, exclusions)
#   differing between runs.
#
# Notes:
#   - Raw RFU data were baseline-corrected and fit using sigmoidal/sliwin models
#   - Amplification efficiency and Cq values were QC-filtered per well
#   - Expression values were normalized using MNE relative to PP2A
#   - Plate-specific scripts are not included to avoid redundancy
#
# Relationship to manuscript analysis:
#   The outputs of this pipeline were merged into `SynthesisMNE.csv`,
#   which is analyzed in `scripts/11b_qPCR.R`.

# import useful libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(chipPCR)
library(qpcR)
library(ggplot2)
library(purrr)


# import dataset
rfu_data <- read_csv("./data/input/qPCR_RFUs.csv", show_col_types = FALSE) %>%
  mutate(Cycle = as.integer(Cycle)) %>%
  mutate(across(-Cycle, as.numeric))


#Change dataset orientation and match plate template
rfu_long <- rfu_data %>%
  pivot_longer(-Cycle, names_to = "Well", values_to = "RFU")

plate_map <- read_csv("./data/input/qPCR_Template.csv", show_col_types = FALSE)

raw <- rfu_long %>%
  left_join(plate_map, by = "Well")


# plot amplification curves
ggplot(raw, aes(x = Cycle, y = RFU, group = Well)) +
  geom_line(linewidth = 0.3, alpha = 0.6, color = "red") +
  theme_bw() +
  labs(title = "Amplification Curves", y = "Fluorescence (RFU)") +
  theme(legend.position = "none")

ggplot(raw, aes(x = Cycle, y = RFU, group = Well)) +
  geom_line(alpha = 0.4, color = "red") +
  facet_wrap(~Target, scales = "free_y") +
  theme_bw() +
  labs(title = "Raw Amplification Curves")


# baseline subtraction
baseline_smooth <- function(cyc, fluo) {
  stopifnot(length(cyc) == length(fluo))
  tibble(Cycle = cyc, RFU = fluo, RFU_corr = fluo - min(fluo, na.rm = TRUE))
}


# per-well fitting function
`%||%` <- function(a, b) if (!is.null(a)) a else b

fit_well <- function(df_one_well) {
  df <- df_one_well %>% arrange(Cycle)
  pp <- baseline_smooth(df$Cycle, df$RFU)
  
  if(max(pp$RFU_corr, na.rm = TRUE) < 1) {
    return(tibble(Cq = NA_real_, Eff = NA_real_, F0 = NA_real_, method = "no_amp"))
  }
  
  # Linear regression of exponential phase
  E_lin <- NA_real_; Cq_lin <- NA_real_; F0_lin <- NA_real_; lin_ok <- FALSE
  try({
    sw <- sliwin(pp$Cycle, log(pmax(pp$RFU_corr, .Machine$double.eps)),
                 width = 4:7, fb = c(2, max(pp$Cycle)-2), simplify = TRUE)
    if(!is.null(sw$par)) {
      E_lin <- as.numeric(sw$par$eff)
      Cq_lin <- as.numeric(sw$par$cp)
      F0_lin <- as.numeric(sw$par$F0)
      lin_ok <- is.finite(E_lin) && E_lin > 1 && E_lin < 2.2 && is.finite(Cq_lin)
    }
  }, silent = TRUE)
  
  # Sigmoidal fit
  fit <- try(pcrfit(data = pp, cyc = 1, fluo = 2, model = l5), silent = TRUE)
  if(inherits(fit, "try-error")) {
    return(tibble(Cq = ifelse(lin_ok, Cq_lin, NA_real_),
                  Eff = ifelse(lin_ok, E_lin, NA_real_),
                  F0 = ifelse(lin_ok, F0_lin, NA_real_),
                  method = ifelse(lin_ok, "sliwin", "failed")))
  }
  
  ef <- try(efficiency(fit, plot = FALSE, type = "Cy0", method = "spline"), silent = TRUE)
  if(inherits(ef, "try-error") || is.null(ef$eff)) {
    ef <- try(efficiency(fit, plot = FALSE, type = "cpD2", method = "spline"), silent = TRUE)
  }
  
  Cq_sig <- suppressWarnings(as.numeric(ef$Cy0 %||% ef$cpD2))
  Eff_sig <- suppressWarnings(as.numeric(ef$eff))
  F0_sig <- suppressWarnings(as.numeric(ef$init2 %||% ef$init1))
  
  if(lin_ok) {
    tibble(Cq = Cq_lin, Eff = E_lin, F0 = F0_lin, method = "sliwin")
  } else {
    tibble(Cq = Cq_sig, Eff = Eff_sig, F0 = F0_sig, method = "sigmoidal")
  }
}


# run per-well modeling
analyze_plate <- function(raw) {
  raw %>%
    group_by(Well) %>%
    arrange(Cycle, .by_group = TRUE) %>%
    summarize(
      data = list(tibble(Cycle = Cycle, RFU = RFU)),
      Sample = first(na.omit(Sample)),
      Target = first(na.omit(Target)),
      Name = first(na.omit(Name)),
      Replicate = first(na.omit(Replicate)),
      HPI = first(na.omit(HPI)),
      Treatment = first(na.omit(Treatment)),
      .groups = "drop"
    ) %>%
    mutate(
      fit = map(data, purrr::possibly(fit_well,
                                      tibble(Cq=NA_real_, Eff=NA_real_, F0=NA_real_, method="failed")))
    ) %>%
    unnest_wider(fit) %>%
    mutate(
      Cq  = map_dbl(Cq, ~ .x[1] %||% NA_real_),
      Eff = map_dbl(Eff, ~ .x[1] %||% NA_real_),
      F0  = map_dbl(F0, ~ .x[1] %||% NA_real_),
      ok_amp = is.finite(Cq) & is.finite(Eff) & Eff > 1.7 & Eff < 2.3,
      note   = ifelse(ok_amp, "ok", "poor fit/efficiency")
    )
}

modeled <- analyze_plate(raw)


# per-well summary of Cq and PCR efficiency
well_summary <- modeled %>%
  dplyr::select(Well, Sample, Target, Name, Replicate, Cq, Eff, ok_amp, note) %>%
  arrange(Sample, Target, Replicate)

write_csv(well_summary,"./data/output/qPCR_Cq_Eff.csv")


#Check number of replicates
replicate_check <- well_summary %>%
  group_by(Name, Target, Sample) %>%
  summarise(
    replicates_present = unique(Replicate),
    n_replicates = n_distinct(Replicate),
    .groups = "drop"
  ) %>%
  mutate(all_three_replicates = (n_replicates == 3 & setequal(replicates_present, 1:3))) %>%
  filter(!all_three_replicates)
replicate_check


#Check Cq homogeneity across technical replicates
techrep_check_all <- well_summary %>%
  group_by(Name, Sample, Target) %>%
  summarise(
    tech1 = first(Cq[Replicate == 1]),
    tech1_well = first(Well[Replicate == 1]),
    tech2 = first(Cq[Replicate == 2]),
    tech2_well = first(Well[Replicate == 2]),
    tech3 = first(Cq[Replicate == 3]),
    tech3_well = first(Well[Replicate == 3]),
    cq_min = min(Cq, na.rm = TRUE),
    cq_max = max(Cq, na.rm = TRUE),
    cq_range = cq_max - cq_min,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    outlier_tech1 = if (!is.na(tech1) & !is.na(tech2) & !is.na(tech3)) (abs(tech1 - tech2) > 0.5 & abs(tech1 - tech3) > 0.5) else FALSE,
    outlier_tech2 = if (!is.na(tech1) & !is.na(tech2) & !is.na(tech3)) (abs(tech2 - tech1) > 0.5 & abs(tech2 - tech3) > 0.5) else FALSE,
    outlier_tech3 = if (!is.na(tech1) & !is.na(tech2) & !is.na(tech3)) (abs(tech3 - tech1) > 0.5 & abs(tech3 - tech2) > 0.5) else FALSE,
    range_flag = cq_range > 0.75,
    note = paste(
      c(
        if(isTRUE(outlier_tech1)) "rep1" else NULL,
        if(isTRUE(outlier_tech2)) "rep2" else NULL,
        if(isTRUE(outlier_tech3)) "rep3" else NULL,
        if(isTRUE(range_flag)) "range" else NULL
      ),
      collapse = ","
    )
  ) %>%
  ungroup() %>%
  dplyr::select(Name, Sample, Target,
                tech1, tech1_well, outlier_tech1,
                tech2, tech2_well, outlier_tech2,
                tech3,tech3_well, outlier_tech3,
                cq_min, cq_max, cq_range, range_flag,
                note) %>%
  arrange(Name, Sample, Target)

write_csv(techrep_check_all,"./data/output/qPCR_checks.csv")


#Eliminate wells based on Maestro analysis + modeling failure + Cq homogeneity checks + blanks
wells_to_remove <- c("F24",
                     "L05","K03",
                     "M17","L06","M19","O18","C04","B06","H22",
                     "A08","L16","M10","B23","C10","J23",
                     "K11","M11","A12","M23","P21","F24","D21",
                     "N22","N23","N24","P22","P23","P24")

well_summary_clean <- well_summary %>%
  mutate(
    Cq = ifelse(Well %in% wells_to_remove, NA, Cq),
    Eff = ifelse(Well %in% wells_to_remove, NA, Eff),
    note = ifelse(Well %in% wells_to_remove, "Excluded", note)
  )

write_csv(well_summary_clean,"./data/output/qPCR_Cq_Eff_cleanup.csv")


#Inspect PCR efficiency distribution
ggplot(well_summary_clean, 
       aes(x = Target, y = Eff)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(title = "PCR Efficiency per Target", y = "Efficiency")


#Adjust wells list for modeled accordingly
modeled_filtered <- modeled %>%
  filter(!Well %in% wells_to_remove)


#Relative quantification - MNE - Normalization against time 0/Mock
calc_MNE <- function(modeled_filtered, baseline_type = c("time0", "mock"), ref_gene = "PP2A") {
  baseline_type <- match.arg(baseline_type)  
  
  #Average NE of reference gene
  ref_NE <- modeled_filtered %>%
    filter(Target == ref_gene) %>%
    mutate(NE_ref = Eff^(-Cq)) %>%
    group_by(Sample, Name) %>%
    summarise(
      NE_ref_mean = exp(mean(log(NE_ref), na.rm = TRUE)),
      .groups = "drop"
    )
  
  #Compute NE for targets
  df <- modeled_filtered %>%
    filter(Target != ref_gene) %>%
    left_join(ref_NE, by = c("Sample", "Name")) %>%
    mutate(NE_target = Eff^(-Cq) / NE_ref_mean)
  
  # Define baseline
  if (baseline_type == "time0") {
    baseline <- df %>%
      filter(HPI == 0) %>%
      group_by(Treatment, Target) %>%
      summarise(NE_baseline = exp(mean(log(NE_target), na.rm = TRUE)), .groups = "drop")
    
    df <- df %>%
      left_join(baseline, by = c("Treatment", "Target"))
    
  } else if (baseline_type == "mock") {
    baseline <- df %>%
      filter(Treatment == "Mock") %>%
      group_by(HPI, Target) %>%
      summarise(NE_baseline = exp(mean(log(NE_target), na.rm = TRUE)), .groups = "drop")
    
    df <- df %>%
      left_join(baseline, by = c("HPI","Target"))
  }
  
  # Fold change relative to baseline
  df <- df %>%
    mutate(NE_fc = NE_target / NE_baseline)
  
  df %>%
    dplyr::select(Well, Sample, Name, Replicate, Target, HPI, Treatment,
                  Cq, Eff, NE_ref_mean, NE_target, NE_baseline, NE_fc)
}

ne_time0 <- calc_MNE(modeled_filtered, baseline_type = "time0")
ne_mock <- calc_MNE(modeled_filtered, baseline_type = "mock")

# summarize and save results
mne_summary_time0 <- ne_time0%>%
  mutate(log_fc = log(NE_fc)) %>%
  group_by(HPI, Treatment, Target, Name, Sample) %>%
  summarise(
    n_replicates = sum(!is.na(NE_fc)),
    wells_used = paste(Well[!is.na(NE_fc)], collapse = ","),
    MNE_fc = exp(mean(log_fc, na.rm = TRUE)),
    se_log = sd(log_fc, na.rm = TRUE) / sqrt(n_replicates),
    se_MNE   = MNE_fc * (exp(se_log) - 1),
    .groups   = "drop"
  )

mne_summary_mock <- ne_mock %>%
  mutate(log_fc = log(NE_fc)) %>%
  group_by(HPI, Treatment, Target, Name, Sample) %>%
  summarise(
    n_replicates = sum(!is.na(NE_fc)),
    wells_used = paste(Well[!is.na(NE_fc)], collapse = ","),
    MNE_fc = exp(mean(log_fc, na.rm = TRUE)),
    se_log = sd(log_fc, na.rm = TRUE) / sqrt(n_replicates),
    se_MNE   = MNE_fc * (exp(se_log) - 1),
    .groups   = "drop"
  )

mne_summary_time0 <- mne_summary_time0 %>%
  mutate(HPI = as.numeric(HPI))
mne_summary_mock <- mne_summary_mock %>%
  mutate(HPI = as.numeric(HPI))

write_csv(ne_time0, "./data/output/qPCR_NE_Time0.csv")
write_csv(ne_mock, "./data/output/qPCR_NE_Mock.csv")
write_csv(mne_summary_time0, "./data/output/qPCR_MNE_summary_Time0.csv")
write_csv(mne_summary_mock, "./data/output/qPCR_MNE_summary_Mock.csv")


# visualize results/generate plots
timepoints <- sort(unique(mne_summary_mock$HPI))
treatments <- unique(mne_summary_time0$Treatment)

plot_per_timepoint <- function(df, timepoint, facet_targets = TRUE) {
  df_sub <- df %>% 
    filter(HPI == timepoint, Target != "PP2A")
  if (nrow(df_sub) == 0) return(NULL)
  
  p <- ggplot(df_sub, aes(x = Treatment, y = MNE_fc, color = Target)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbar(
      aes(ymin = MNE_fc - se_MNE, ymax = MNE_fc + se_MNE),
      width = 0.2,
      position = position_dodge(width = 0.5)
    ) +
    theme_bw() +
    labs(
      title = paste("Expression at", timepoint, "hpi"),
      y = "Relative Expression (MNE)",
      x = "Treatment"
    )
  
  if (facet_targets) {
    p <- p + facet_wrap(~Target, scales = "free_y")
  }
  
  p
}

plot_per_treatment <- function(df, treatment, facet_targets = TRUE) {
  df_sub <- df %>% 
    filter(Treatment == treatment, Target != "PP2A")
  if (nrow(df_sub) == 0) return(NULL)
  
  p <- ggplot(df_sub, aes(x = HPI, y = MNE_fc, color = Target, group = Target)) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = MNE_fc - se_MNE, ymax = MNE_fc + se_MNE),
      width = 0.2
    ) +
    geom_line() +
    theme_bw() +
    labs(
      title = paste("Expression over time - ", treatment),
      y = "Relative Expression (MNE)",
      x = "Time (hpi)"
    )
  
  if (facet_targets) {
    p <- p + facet_wrap(~Target, scales = "free_y")
  }
  
  p
}

plots_timepoints <- lapply(timepoints, function(tp) plot_per_timepoint(mne_summary_mock, tp))
plots_treatments <- lapply(treatments, function(tr) plot_per_treatment(mne_summary_time0, tr))
for (p in plots_timepoints) if (!is.null(p)) print(p)
for (p in plots_treatments) if (!is.null(p)) print(p)


# save as PDF with multiple pages
save_path <- "./figures/qPCR_plots_mne_time0_treatment.pdf"
pdf(save_path, width = 12, height = 8)
for (p in plots_treatments) {
  if (!is.null(p)) print(p)
}
dev.off()

save_path <- "./figures/qPCR_plots_mne_mock_timepoint.pdf"
pdf(save_path, width = 12, height = 8)
for (p in plots_timepoints) {
  if (!is.null(p)) print(p)
}
dev.off()