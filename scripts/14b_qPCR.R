# 11_qPCR.R
#
# Purpose:
#   Analyze qPCR expression data from infected Arabidopsis seedlings using
#   MNE-normalized fold change (MNE_fc) values. The script summarizes expression
#   across biological replicates, generates time-course and treatment-specific
#   visualizations, and performs adaptive statistical testing for each immune
#   marker gene.
#
# Input:
#   - SynthesisMNE.csv
#     Raw qPCR output containing MNE_fc values, targets, treatments, timepoints,
#     and biological replicate information.
#
# Outputs:
#   - SynthesisMNE_summary.csv
#       Summary table of mean expression and standard error per target, treatment,
#       and timepoint.
#   - Multi-page PDF figures showing:
#       • Expression by treatment at each timepoint
#       • Expression over time per treatment
#       • Paper-ready plots for selected immune marker genes
#
# Methods:
#   - Log-transformation of MNE_fc values prior to averaging
#   - Back-transformation to geometric means and propagated SE
#   - Adaptive statistical testing per target and timepoint:
#       • ANOVA + Tukey HSD
#       • Welch ANOVA + Games–Howell
#       • Kruskal–Wallis + Dunn’s test
#     Selection based on Shapiro–Wilk normality and Levene’s variance tests
#
# Dependencies:
#   tidyverse, chipPCR, qpcR, ggplot2, purrr, broom,
#   FSA, car, rstatix, dunn.test, viridis
#


#Import useful libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(chipPCR)
library(qpcR)
library(ggplot2)
library(purrr)
library(broom)
library(FSA)
library(carData)
library(car)
library(rstatix)
library(dunn.test)
library(stringr)
library(viridisLite)
library(viridis)


# import dataset
Synthesis <- read_csv("./data/input/qPCR_SynthesisMNE.csv", show_col_types = FALSE)

# calculate mean MNE_fc + SE
Synthesis_summary <- Synthesis %>%
  mutate(log_MNE_fc = log(MNE_fc)) %>%
  group_by(HPI, Treatment, Name, Target, Baseline) %>%
  summarise(
    n_replicates = sum(!is.na(Biorep)),
    Samples = paste(Sample[!is.na(MNE_fc)], collapse = ","),
    mean_log_MNE_fc = mean(log_MNE_fc, na.rm = TRUE),
    MNE_fc_global = exp(mean_log_MNE_fc),
    se_log_MNE_fc = sd(log_MNE_fc, na.rm = TRUE)/sqrt(n_replicates),
    se_MNE_fc_global = MNE_fc_global * (exp(se_log_MNE_fc) - 1),
    .groups   = "drop"
  )

write_csv(Synthesis_summary,"./data/output/qPCR_SynthesisMNE_summary.csv")


# visualize results/generate plots - All data
Synthesis_summary_mock <- subset(Synthesis_summary, Synthesis_summary$Baseline == "Mock")
timepoints <- sort(unique(Synthesis_summary_mock$HPI))

Synthesis_summary_time0 <- subset(Synthesis_summary, Synthesis_summary$Baseline == "Time0")
treatments <- unique(Synthesis_summary_time0$Treatment)

plot_per_timepoint <- function(df, timepoint, facet_targets = TRUE) {
  df_sub <- df %>% 
    filter(HPI == timepoint)
  if (nrow(df_sub) == 0) return(NULL)
  
  p <- ggplot(df_sub, aes(x = Treatment, y = MNE_fc_global, color = Target)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbar(
      aes(ymin = MNE_fc_global - se_MNE_fc_global, ymax = MNE_fc_global + se_MNE_fc_global),
      width = 0.2,
      position = position_dodge(width = 0.5)
    ) +
    theme_bw() +
    labs(
      title = paste("Expression at", timepoint, "hpi"),
      y = "Relative Expression (MNE_fc)",
      x = "Treatment"
    )
  
  if (facet_targets) {
    p <- p + facet_wrap(~Target, scales = "free_y")
  }
  
  p
}

plot_per_treatment <- function(df, treatment, facet_targets = TRUE) {
  df_sub <- df %>% 
    filter(Treatment == treatment)
  if (nrow(df_sub) == 0) return(NULL)
  
  p <- ggplot(df_sub, aes(x = HPI, y = MNE_fc_global, color = Target, group = Target)) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = MNE_fc_global - se_MNE_fc_global, ymax = MNE_fc_global + se_MNE_fc_global),
      width = 0.2
    ) +
    geom_line() +
    theme_bw() +
    labs(
      title = paste("Expression over time", treatment),
      y = "Relative Expression (MNE_fc)",
      x = "Time (hpi)"
    )
  
  if (facet_targets) {
    p <- p + facet_wrap(~Target, scales = "free_y")
  }
  
  p
}

plots_timepoints <- lapply(timepoints, function(tp) plot_per_timepoint(Synthesis_summary_mock, tp))
plots_treatments <- lapply(treatments, function(tr) plot_per_treatment(Synthesis_summary_time0, tr))
for (p in plots_timepoints) if (!is.null(p)) print(p)
for (p in plots_treatments) if (!is.null(p)) print(p)


# save as PDF
pdf("./figures/qPCR_time0_treatment_all_data.pdf", width = 12, height = 8)
for (p in plots_treatments) {
  if (!is.null(p)) print(p)
}
dev.off()

pdf("./figures/qPCR_mock_timepoints_all_data.pdf", width = 12, height = 8)
for (p in plots_timepoints) {
  if (!is.null(p)) print(p)
}
dev.off()


# select useful data & create subsets x target
MNE <- subset(Synthesis, Synthesis$Baseline=="Mock" & Synthesis$Treatment!= "Mock" & Synthesis$Treatment!= "OBC")
MNE_summary <- subset(Synthesis_summary, Synthesis_summary$Baseline=="Mock" & Synthesis_summary$Treatment!= "Mock" & Synthesis_summary$Treatment!= "OBC")

FRK1 <- subset(MNE_summary, MNE_summary$Target== "FRK1" & MNE_summary$HPI!= "72" & MNE_summary$HPI!= "48" & MNE_summary$HPI!= "24")
WRKY29 <- subset(MNE_summary, MNE_summary$Target== "WRKY29" & MNE_summary$HPI!= "72" & MNE_summary$HPI!= "48" & MNE_summary$HPI!= "24")
NHL10 <- subset(MNE_summary, MNE_summary$Target== "NHL10")
PAD3 <- subset(MNE_summary, MNE_summary$Target== "PAD3")
PDF1.2 <- subset(MNE_summary, MNE_summary$Target== "PDF1.2" & MNE_summary$HPI!= "2" & MNE_summary$HPI!= "6")
PR1 <- subset(MNE_summary, MNE_summary$Target== "PR1" & MNE_summary$HPI!= "2" & MNE_summary$HPI!= "6")


# visualize results/generate plots - x paper
plotFRK1 <- ggplot(FRK1, aes(x = HPI, y = MNE_fc_global, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MNE_fc_global - se_MNE_fc_global,
        ymax = MNE_fc_global + se_MNE_fc_global),
    width = 0.2
  ) +
  theme_bw() + scale_color_viridis_d(option = "D") +
  labs(
    title = "FRK1 baseline Mock",
    y = "Relative Expression (MNE_fc) +/- SE",
    x = "HPI"
  )
plotFRK1

plotWRKY29 <- ggplot(WRKY29, aes(x = HPI, y = MNE_fc_global, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MNE_fc_global - se_MNE_fc_global,
        ymax = MNE_fc_global + se_MNE_fc_global),
    width = 0.2
  ) +
  theme_bw() + scale_color_viridis_d(option = "D") +
  labs(
    title = "WRKY29 baseline Mock",
    y = "Relative Expression (MNE_fc) +/- SE",
    x = "HPI"
  )
plotWRKY29

plotNHL10 <- ggplot(NHL10, aes(x = HPI, y = MNE_fc_global, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MNE_fc_global - se_MNE_fc_global,
        ymax = MNE_fc_global + se_MNE_fc_global),
    width = 0.2
  ) +
  theme_bw() + scale_color_viridis_d(option = "D") +
  labs(
    title = "NHL10 baseline Mock",
    y = "Relative Expression (MNE_fc) +/- SE",
    x = "HPI"
  )
plotNHL10

plotPAD3 <- ggplot(PAD3, aes(x = HPI, y = MNE_fc_global, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MNE_fc_global - se_MNE_fc_global,
        ymax = MNE_fc_global + se_MNE_fc_global),
    width = 0.2
  ) +
  theme_bw() + scale_color_viridis_d(option = "D") +
  labs(
    title = "PAD3 baseline Mock",
    y = "Relative Expression (MNE_fc) +/- SE",
    x = "HPI"
  )
plotPAD3

plotPDF1.2 <- ggplot(PDF1.2, aes(x = HPI, y = MNE_fc_global, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MNE_fc_global - se_MNE_fc_global,
        ymax = MNE_fc_global + se_MNE_fc_global),
    width = 0.2
  ) +
  theme_bw() + scale_color_viridis_d(option = "D") +
  labs(
    title = "PDF1.2 baseline Mock",
    y = "Relative Expression (MNE_fc) +/- SE",
    x = "HPI"
  )
plotPDF1.2

plotPR1 <- ggplot(PR1, aes(x = HPI, y = MNE_fc_global, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = MNE_fc_global - se_MNE_fc_global,
        ymax = MNE_fc_global + se_MNE_fc_global),
    width = 0.2
  ) +
  theme_bw() + scale_color_viridis_d(option = "D") +
  labs(
    title = "PR1 baseline Mock",
    y = "Relative Expression (MNE_fc) +/- SE",
    x = "HPI"
  )
plotPR1


# Save as PDF with multiple pages
plots <- list(
  plotFRK1,
  plotWRKY29,
  plotNHL10,
  plotPAD3,
  plotPDF1.2,
  plotPR1
)
pdf("./figures/qPCR_targets_x_paper.pdf", width = 12, height = 8)
for (p in plots) {
  if (!is.null(p)) print(p)
}
dev.off()


# create new subsets x target with all bioreps for stats computation
FRK1 <- subset(MNE, MNE$Target== "FRK1" & MNE$HPI!= "72" & MNE$HPI!= "48" & MNE$HPI!= "24")
WRKY29 <- subset(MNE, MNE$Target== "WRKY29" & MNE$HPI!= "72" & MNE$HPI!= "48" & MNE$HPI!= "24")
NHL10 <- subset(MNE, MNE$Target== "NHL10")
PAD3 <- subset(MNE, MNE$Target== "PAD3")
PDF1.2 <- subset(MNE, MNE$Target== "PDF1.2" & MNE$HPI!= "2" & MNE$HPI!= "6")
PR1 <- subset(MNE, MNE$Target== "PR1" & MNE$HPI!= "2" & MNE$HPI!= "6")


# stats
run_adaptive_stats <- function(df, target_name) {
  df <- df %>%
    filter(!is.na(MNE_fc)) %>%
    mutate(Treatment = factor(Treatment))
  
  mod <- lm(MNE_fc ~ Treatment, data = df)
  resid_vals <- residuals(mod)
  
  sh_p <- if(sd(resid_vals) == 0) NA else shapiro.test(resid_vals)$p.value
  lev_p <- if(sd(df$MNE_fc) == 0) NA else leveneTest(MNE_fc ~ Treatment, data = df)$`Pr(>F)`[1]
  
  if(!is.na(sh_p) && sh_p > 0.05 && !is.na(lev_p) && lev_p > 0.05){
    aov_res <- aov(MNE_fc ~ Treatment, data = df)
    global_p <- summary(aov_res)[[1]]$`Pr(>F)`[1]
    test_used <- "ANOVA"
    posthoc <- if(!is.na(global_p) && global_p < 0.05){
      broom::tidy(TukeyHSD(aov_res)) %>%
        select(term, group1, group2, p.adj)
    } else NULL
    
  } else if(!is.na(sh_p) && sh_p > 0.05 && !is.na(lev_p) && lev_p <= 0.05){
    welch_res <- oneway.test(MNE_fc ~ Treatment, data = df, var.equal = FALSE)
    global_p <- welch_res$p.value
    test_used <- "Welch ANOVA"
    posthoc <- if(!is.na(global_p) && global_p < 0.05){
      gh <- games_howell_test(df, MNE_fc ~ Treatment)
      if("group1" %in% colnames(gh)){
        gh %>% select(group1, group2, p.adj)
      } else if("comparison" %in% colnames(gh)){
        gh %>% tidyr::separate(comparison, into = c("group1","group2"), sep = " - ") %>% select(group1, group2, p.adj)
      } else gh
    } else NULL
    
  } else {
    kw_res <- kruskal_test(MNE_fc ~ Treatment, data = df)
    global_p <- kw_res$p
    test_used <- "Kruskal-Wallis"
    posthoc <- if(!is.na(global_p) && global_p < 0.05){
      dunnTest(MNE_fc ~ Treatment, data = df, method = "bonferroni")$res %>%
        rename(p.adj = P.adj) %>%
        tidyr::separate(Comparison, into = c("group1","group2"), sep = " - ") %>%
        select(group1, group2, p.adj)
    } else NULL
  }
  
  tibble(
    Target = target_name,
    HPI = unique(df$HPI),
    Test = test_used,
    Global_p = global_p,
    Posthoc = list(posthoc)
  )
}

analyze_target <- function(target_df, target_name){
  target_df %>%
    filter(!is.na(HPI)) %>%
    group_by(HPI) %>%
    group_modify(~ run_adaptive_stats(.x, target_name)) %>%
    ungroup()
}

FRK1_stats   <- analyze_target(FRK1, "FRK1")
WRKY29_stats <- analyze_target(WRKY29, "WRKY29")
NHL10_stats  <- analyze_target(NHL10, "NHL10")
PAD3_stats   <- analyze_target(PAD3, "PAD3")
PDF1.2_stats <- analyze_target(PDF1.2, "PDF1.2")
PR1_stats    <- analyze_target(PR1, "PR1")

FRK1_stats
WRKY29_stats
NHL10_stats
PAD3_stats
PDF1.2_stats
PR1_stats