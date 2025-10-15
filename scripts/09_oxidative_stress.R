# 09_oxidative_stress.R
#
# Purpose:
#   Test whether O-antigen (OPS) mutants differ from WT in oxidative stress tolerance.
#   Computes per-well AUCs, normalizes to 0 mM, estimates per-plate WT MIC,
#   scales doses by WT MIC, fits mixed models, extracts emmeans (AUC on response scale),
#   and saves MICs, WT-vs-mutant contrasts, and the final figure.
#
# Inputs  (place files here):
#   data/09_oxidative_stress/input/
#     - OD600_R12.csv ... OD600_R16.csv
#     - layout_R12.csv ... layout_R16.csv
#
# Outputs:
#   figures/09_oxidative_stress.pdf
#   data/09_oxidative_stress/output/WT_MIC_by_replicate.csv
#   data/09_oxidative_stress/output/emmeans_contrasts_vsWT_scaledMIC.csv
#
# Dependencies:
#   ggplot2, viridis, dplyr, tidyr, readr, stringr, forcats, janitor, purrr,
#   lme4, lmerTest, emmeans
#
# Notes:
#   - Dose is scaled per plate by the WT MIC estimated from AUC_rel0 vs H2O2.
#   - Final outputs match the Results text (equal MICs; no significant AUC differences).

library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(janitor)
library(purrr)
library(lme4)
library(lmerTest)
library(emmeans)

# ------------------ paths ------------------
input_dir <- "data/09_oxidative_stress/input"
out_dir   <- "data/09_oxidative_stress/output"
out_fig   <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

# Optional QC saves (raw curves/AUC per plate). Set TRUE only if you want extras.
SAVE_QC <- FALSE
qc_dir  <- "figures/09_oxidative_stress_QC"
if (SAVE_QC) dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------ helpers ------------------
normalize_well <- function(x) {
  x <- toupper(trimws(x))
  gsub("^([A-H])0?([1-9]|1[0-2])$", "\\1\\2", x)
}

trapezoid_auc <- function(t, y) {
  o <- order(t); t <- t[o]; y <- y[o]
  if (length(t) < 2) return(NA_real_)
  sum(diff(t) * (head(y, -1) + tail(y, 1)) / 2, na.rm = TRUE)
}

# Estimate MIC by linear interpolation of WT AUC_rel0 vs conc (threshold 0.1)
estimate_mic_interp <- function(df_one_rep) {
  df <- df_one_rep %>%
    group_by(H2O2conc) %>%
    summarise(y = mean(AUC_rel0, na.rm = TRUE), .groups = "drop") %>%
    arrange(H2O2conc)
  if (nrow(df) < 2 || all(df$y > 0.1)) return(NA_real_)
  idx <- which(df$y[-1] <= 0.1 & df$y[-length(df$y)] > 0.1)
  if (length(idx) == 0) return(min(df$H2O2conc[df$y <= 0.1], na.rm = TRUE))
  i <- idx[1]
  x0 <- df$H2O2conc[i]; x1 <- df$H2O2conc[i + 1]
  y0 <- df$y[i];        y1 <- df$y[i + 1]
  x0 + (0.1 - y0) * (x1 - x0) / (y1 - y0)
}

read_layout <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  nm  <- names(raw)
  
  # Matrix-style layout (row + columns 1..12)
  if (any(grepl("^pcr\\.?plate$", nm, ignore.case = TRUE)) && all(as.character(1:12) %in% nm)) {
    out <- raw %>%
      rename(row = 1) %>%
      pivot_longer(-row, names_to = "col", values_to = "entry") %>%
      mutate(
        well   = normalize_well(paste0(row, gsub("\\D", "", col))),
        entry  = as.character(entry),
        is_blank = grepl("^plainLB$", entry, ignore.case = TRUE),
        strain   = if_else(is_blank, NA_character_, sub("\\s*\\+.*$", "", entry)),
        conc_str = if_else(is_blank, NA_character_, sub("^.*\\+\\s*", "", entry)),
        H2O2conc = suppressWarnings(as.numeric(conc_str)),
        strain   = case_when(
          is.na(strain) ~ NA_character_,
          grepl("^p25C2$", strain, ignore.case = TRUE) ~ "WT",
          TRUE ~ gsub("^∆", "Δ", strain)
        )
      ) %>%
      select(well, strain, H2O2conc) %>%
      filter(!is.na(strain))
    return(out)
  }
  
  # Tidy fallback: first three cols well/strain/H2O2conc
  raw %>%
    clean_names() %>%
    rename(well = 1, strain = 2, H2O2conc = 3) %>%
    mutate(
      well   = normalize_well(well),
      strain = if_else(grepl("^p25C2$", strain, ignore.case = TRUE), "WT",
                       gsub("^∆", "Δ", as.character(strain))),
      H2O2conc = as.numeric(H2O2conc)
    ) %>%
    filter(!is.na(strain))
}

read_OD <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  tcol  <- names(df)[grepl("^time", names(df), ignore.case = TRUE)][1]
  if (is.na(tcol)) stop("Time column not found in: ", path)
  wells <- names(df)[grepl("^[A-H](?:0?[1-9]|1[0-2])$", names(df))]
  if (length(wells) == 0) stop("No well columns detected in: ", path)
  
  df %>%
    rename(time_s = all_of(tcol)) %>%
    mutate(time_h = as.numeric(time_s) / 3600) %>%
    pivot_longer(all_of(wells), names_to = "well", values_to = "OD600") %>%
    mutate(
      well  = normalize_well(well),
      time  = time_h,
      OD600 = as.numeric(OD600)
    ) %>%
    select(time, well, OD600)
}

# ------------------ ingest all replicates ------------------
rep_ids <- c("R12", "R13", "R15", "R16")

long_all <- map_dfr(rep_ids, function(r) {
  layout_path <- file.path(input_dir, paste0("layout_", r, ".csv"))
  od_path     <- file.path(input_dir, paste0("OD600_", r, ".csv"))
  if (!file.exists(layout_path) || !file.exists(od_path)) {
    warning(sprintf("Skipping %s: missing files.", r))
    return(tibble())
  }
  L <- read_layout(layout_path) %>%
    mutate(well = normalize_well(well))
  O <- read_OD(od_path)
  inner_join(O, L, by = "well") %>%
    mutate(replicate = r)
})

stopifnot(nrow(long_all) > 0)

# ------------------ per-well AUC (0–24 h) + normalize to 0 mM ------------------
auc_per_well <- long_all %>%
  filter(time >= 0, time <= 24) %>%
  group_by(replicate, well, strain, H2O2conc) %>%
  summarise(AUC = trapezoid_auc(time, OD600), .groups = "drop")

ref0 <- auc_per_well %>%
  filter(H2O2conc == 0) %>%
  group_by(replicate, strain) %>%
  summarise(AUC0 = mean(AUC, na.rm = TRUE), .groups = "drop")

auc_norm <- auc_per_well %>%
  left_join(ref0, by = c("replicate", "strain")) %>%
  mutate(AUC_rel0 = AUC / AUC0) %>%
  filter(is.finite(AUC_rel0))

# ------------------ estimate WT MIC per replicate + scale dose ------------------
wt_mic <- auc_norm %>%
  filter(strain == "WT") %>%
  group_by(replicate, H2O2conc) %>%
  summarise(AUC_rel0 = mean(AUC_rel0, na.rm = TRUE), .groups = "drop") %>%
  group_by(replicate) %>%
  group_modify(~tibble(MIC_WT = estimate_mic_interp(.x))) %>%
  ungroup()

# Fallback: if threshold not crossed, use max tested conc
fallback <- auc_norm %>% filter(strain == "WT") %>%
  group_by(replicate) %>% summarise(maxC = max(H2O2conc, na.rm = TRUE), .groups = "drop")

wt_mic <- wt_mic %>%
  left_join(fallback, by = "replicate") %>%
  mutate(MIC_WT = ifelse(is.na(MIC_WT), maxC, MIC_WT)) %>%
  select(-maxC)

# Save MICs (reported in Results)
write_csv(wt_mic, file.path(out_dir, "WT_MIC_by_replicate.csv"))

auc_scaled <- auc_norm %>%
  left_join(wt_mic, by = "replicate") %>%
  mutate(H2O2_scaled = H2O2conc / MIC_WT) %>%
  filter(is.finite(H2O2_scaled))

# ------------------ mixed model on scaled dose + emmeans ------------------
model_df <- auc_scaled %>%
  mutate(
    strain = fct_relevel(factor(strain), "WT"),
    logAUC = log(AUC_rel0 + 1e-6)
  )

fit <- lmer(logAUC ~ strain * H2O2_scaled + (1 | replicate), data = model_df)

# emmeans at standard scaled doses (response scale)
dose_grid <- c(0, 0.25, 0.5, 0.75, 1.0, 1.25)
emm <- emmeans(fit, ~ strain | H2O2_scaled,
               at = list(H2O2_scaled = dose_grid),
               type = "response")

# contrasts vs WT (+ FDR)
contr_df <- as.data.frame(summary(contrast(emm, method = "trt.vs.ctrl", ref = "WT"),
                                  type = "response"))
val_col <- if ("response" %in% names(contr_df)) "response" else if ("ratio" %in% names(contr_df)) "ratio" else "estimate"
LCL_col <- if ("lower.CL" %in% names(contr_df)) "lower.CL" else if ("LCL" %in% names(contr_df)) "LCL" else NA_character_
UCL_col <- if ("upper.CL" %in% names(contr_df)) "upper.CL" else if ("UCL" %in% names(contr_df)) "UCL" else NA_character_

contr_out <- contr_df %>%
  mutate(
    ratio    = .data[[val_col]],
    lower.CL = if (!is.na(LCL_col)) .data[[LCL_col]] else NA_real_,
    upper.CL = if (!is.na(UCL_col)) .data[[UCL_col]] else NA_real_,
    p.adj    = p.adjust(p.value, method = "fdr")
  ) %>%
  select(H2O2_scaled, contrast, ratio, lower.CL, upper.CL, p.value, p.adj)

# Save contrasts (reported in Results)
write_csv(contr_out, file.path(out_dir, "emmeans_contrasts_vsWT_scaledMIC.csv"))

# ------------------ final figure ------------------
emm_means <- as.data.frame(emm)
mean_col  <- if ("response" %in% names(emm_means)) "response" else "emmean"
LCLm_col  <- if ("lower.CL" %in% names(emm_means)) "lower.CL" else "LCL"
UCLm_col  <- if ("upper.CL" %in% names(emm_means)) "upper.CL" else "UCL"

emm_plot <- emm_means %>%
  mutate(
    mean_bt  = .data[[mean_col]],
    lower_bt = if (LCLm_col %in% names(.)) .data[[LCLm_col]] else NA_real_,
    upper_bt = if (UCLm_col %in% names(.)) .data[[UCLm_col]] else NA_real_
  )

# order strains for plotting
emm_plot$strain <- factor(emm_plot$strain,
  levels = c("epsE", "rmlC", "spsA", "tagG", "tagH", "wfgD", "WT"))
# plot
p <- ggplot(emm_plot, aes(x = H2O2_scaled, y = mean_bt, color = strain, group = strain)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_bt, ymax = upper_bt), width = 0.03) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 1, linetype = 3) +
  scale_color_viridis_d(option = "D") +
  theme_bw(base_size = 16) +
  labs(
    x = "Scaled dose (H2O2 / WT MIC)",
    y = "AUC (rel. 0 mM; emmeans, response scale)",
    title = "Strain responses after per-plate WT MIC scaling")
p

ggsave(file.path(out_fig, "09_oxidative_stress.pdf"), p, width = 9, height = 6, dpi = 300)
