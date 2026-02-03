# 14_AMP.R
#
# Purpose:
#   Analyze antimicrobial peptide (AMP) disk diffusion assays for wild-type and
#   O-antigen mutant strains. The script summarizes inhibition measurements across
#   replicates, generates dose–response plots by compound, and performs adaptive
#   statistical tests comparing strains at each dose.
#
# Input:
#   - AMP_Synthesis.csv
#     Raw disk diffusion dataset containing inhibition measurements (MEANAREA),
#     strain identities, compound identity, disk number, and applied quantity.
#
# Outputs:
#   - Synthesis_analysis.csv
#       Summary table of mean inhibition area (radius) ± standard error for each
#       strain × compound × dose combination.
#   - Overall_analysis.pdf
#       Multi-page PDF containing dose–response plots for each compound (PMB, PME).
#   - AMP_PMB_stats.csv
#       Per-dose global test results comparing strains for PMB (adaptive testing).
#   - AMP_PME_stats.csv
#       Per-dose global test results comparing strains for PME (adaptive testing).
#
# Methods:
#   - Replicate summarization: mean(MEANAREA) and SE = sd/sqrt(n) by STRAIN × COMPOUND × DISK NUMBER × DOSE
#   - Visualization: line plots of inhibition radius vs dose per compound, colored by strain
#   - Adaptive statistical testing at each dose:
#       • ANOVA + Tukey HSD (if residuals normal and variances equal)
#       • Welch ANOVA + Games–Howell (if residuals normal but variances unequal)
#       • Kruskal–Wallis + Dunn’s test (if residuals non-normal)
#     Selection based on Shapiro–Wilk normality and Levene’s variance tests.
#
# Dependencies:
#   tidyverse, dplyr, tidyr, ggplot2, purrr, car, rstatix, broom, viridis
#
# Notes:
#   - Quantities are coerced to numeric (`QUANTITY (ug)`).
#   - For per-dose tests, doses with zero variance across replicates are removed
#     (sd(MEANAREA) == 0) to avoid invalid model fits.

#Import useful libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(chipPCR)
library(qpcR)
library(ggplot2)
library(purrr)
library(car)
library(rstatix)
library(broom)
library(viridisLite)
library(viridis)

# import dataset
Dataset <- read_csv("./data/input/AMP_synthesis.csv", show_col_types = FALSE)

# calculate mean inhibition area + SE
Synthesis <- Dataset %>%
  group_by(STRAIN, COMPOUND, `DISK NUMBER`, `QUANTITY (ug)`) %>%
  summarise(
    N_REPLICATES = sum(!is.na(MEANAREA)),
    MEAN_AREA = mean(MEANAREA, na.rm = TRUE),
    SE_AREA = sd(MEANAREA, na.rm = TRUE)/sqrt(N_REPLICATES),
    .groups   = "drop"
  )
Synthesis$`QUANTITY (ug)` <- as.numeric(Synthesis$`QUANTITY (ug)`)

write_csv(Synthesis, "./data/output/AMP_analysis.csv")

# visualize results/generate plots
plot_line_compound <- function(df, compound_name) {
  df %>%
    filter(COMPOUND == compound_name) %>%
    ggplot(aes(x = `QUANTITY (ug)`, y = MEAN_AREA,
               color = STRAIN, group = STRAIN)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(aes(ymin = MEAN_AREA - SE_AREA,
                      ymax = MEAN_AREA + SE_AREA),
                  width = 0.05, linewidth = 0.6) +
    theme_bw(base_size = 14) +
    scale_color_viridis_d(option="D") +
    labs(
      x = "Quantity (µg)",
      y = "Mean inhibition radius (mm) ± SE",
      color = "Strain",
      title = paste("Inhibition radius -", compound_name)
    ) +
    theme(panel.grid.minor = element_blank())
}

plot_PMB <- plot_line_compound(Synthesis, "PMB")
plot_PME <- plot_line_compound(Synthesis, "PME")

# save as PDF 
pdf("./figures/AMP.pdf", width = 12, height = 8)
print(plot_PMB)
print(plot_PME)
dev.off()

# stats PMB/PME
run_adaptive_test <- function(df) {
  mod <- lm(MEANAREA ~ STRAIN, data = df)
  resid_vals <- residuals(mod)
  
  sh_p <- if(sd(resid_vals) == 0) NA else shapiro.test(resid_vals)$p.value
  lev_p <- if(sd(df$MEANAREA) == 0) NA else leveneTest(MEANAREA ~ STRAIN, data = df)$`Pr(>F)`[1]
  
  if(!is.na(sh_p) && sh_p > 0.05 && !is.na(lev_p) && lev_p > 0.05) {
    aov_res <- aov(MEANAREA ~ STRAIN, data = df)
    global_p <- summary(aov_res)[[1]]$`Pr(>F)`[1]
    test_used <- "ANOVA"
    posthoc <- if(global_p < 0.05) TukeyHSD(aov_res) %>% broom::tidy() else NULL
    
  } else if(!is.na(sh_p) && sh_p > 0.05 && !is.na(lev_p) && lev_p <= 0.05) {
    welch_res <- oneway_test(MEANAREA ~ STRAIN, data = df, var.equal = FALSE)
    global_p <- welch_res$p
    test_used <- "Welch ANOVA"
    posthoc <- if(global_p < 0.05) games_howell_test(df, MEANAREA ~ STRAIN) else NULL
    
  } else {
    kw_res <- kruskal_test(MEANAREA ~ STRAIN, data = df)
    global_p <- kw_res$p
    test_used <- "Kruskal-Wallis"
    posthoc <- if(global_p < 0.05) dunn_test(df, MEANAREA ~ STRAIN, p.adjust.method = "bonferroni") else NULL
  }
  
  tibble(
    QUANTITY = unique(df$`QUANTITY (ug)`),
    Test = test_used,
    Global_p = global_p,
    Posthoc = list(posthoc)
  )
}

analyze_compound <- function(compound_name) {
  data_comp <- Dataset %>%
    filter(COMPOUND == compound_name) %>%
    mutate(`QUANTITY (ug)` = as.numeric(`QUANTITY (ug)`)) %>%
    group_by(`QUANTITY (ug)`) %>%
    filter(sd(MEANAREA, na.rm = TRUE) > 0) %>%
    ungroup()
  
  results <- data_comp %>%
    group_by(`QUANTITY (ug)`) %>%
    group_modify(~ run_adaptive_test(.x)) %>%
    ungroup()
  
  results
}

PMB_stats <- analyze_compound("PMB")
PME_stats <- analyze_compound("PME")

PMB_stats <- PMB_stats %>% rename(QUANTITY = `QUANTITY (ug)`)
PME_stats <- PME_stats %>% rename(QUANTITY = `QUANTITY (ug)`)

write_csv(PMB_stats, "./data/output/AMP_PMB_stats.csv")
write_csv(PME_stats, "./data/output/AMP_PME_stats.csv")
