# ============================================================
# 12_invitro_growth.R
#
# Purpose:
#   Analyze in vitro growth dynamics of wild-type and O-antigen mutants 
#   grown overnight at an initial OD600 of 0.01. 
#   - Fit logistic growth models using Growthcurver
#   - Extract per-strain growth parameters (k, r, AUC, etc.)
#   - Compare mutants vs WT using appropriate tests (ANOVA/Tukey, 
#     Welch/Games-Howell, or Kruskal-Wallis/Dunn based on assumptions)
#   - Output summary statistics and growth parameter tables
#
# Inputs:
#   data/12_invitro_growth/Growth_curves_synthesis.csv
#     Required columns:
#       - TIME          (numeric, hours)
#       - ABS           (absorbance at 600 nm)
#       - STRAIN        (strain name)
#       - GENOTYPE      (mutation identifier)
#       - METHOD        (e.g., “Overnight”)
#       - OD            (initial optical density)
#
# Outputs:
#   figures/12_growthcurves.pdf
#   data/12_invitro_growth/output/growthcurveR_data.csv
#   data/12_invitro_growth/output/growthcurve_stats.csv
#
# Dependencies:
#   tidyverse, ggplot2, viridis, growthcurver, purrr, broom, car, rstatix, dunn.test
#
# Notes:
#   - Fits logistic models with SummarizeGrowth() per strain replicate.
#   - Automatically selects post hoc tests based on normality (Shapiro) 
#     and homogeneity (Levene) results.
#   - Significance codes follow conventional thresholds (* p<0.05, ** p<0.01, *** p<0.001).
# ============================================================

library(tidyverse)
library(ggplot2)
library(viridis)
library(growthcurver)
library(purrr)
library(broom)
library(car)
library(rstatix)
library(dunn.test)


# Import and organize data
dSynthesis <- read.csv("./data/12_invitro_growth/Growth_curves_synthesis.csv")
dSynthesis$TIME <- as.numeric(dSynthesis$TIME)
dSynthesisON <- subset(dSynthesis,dSynthesis$METHOD=="Overnight")
dSynthesisON0.01 <- subset(dSynthesisON,dSynthesisON$OD=="0.01")
SINGLE <- subset(dSynthesisON0.01,dSynthesisON0.01$GENOTYPE!="Double")
PAPER <- subset(SINGLE,SINGLE$STRAIN!="OBC" & SINGLE$STRAIN!="rmlC-FS")

# Plot all strains - optical density
ggplot(PAPER, aes(TIME, ABS, color=STRAIN)) + 
  geom_point() + 
  theme_bw() + scale_color_viridis_d(option = "D") + xlab("Time (hours)") + ylab("Absorbance (600nm)") + ggtitle("Single mutants")

SMABS <- ggplot(PAPER, aes(TIME, ABS, color=STRAIN)) + 
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(fun.data = mean_se, geom = "ribbon", aes(fill = STRAIN), alpha = 0.2, color = NA) + 
  theme_bw() + scale_color_viridis_d(option = "D") + scale_fill_viridis_d(option = "D") + xlab("Time (hours)") + ylab("Absorbance (600nm)") + ggtitle("Single mutants")

ggsave("./figures/12_growthcurves.pdf", SMABS, width = 10, height = 6, dpi = 300)


# Fit growth curves individually - optical density
Abs_fits <- PAPER %>%
  group_by(DATE, STRAIN, GENOTYPE) %>%
  arrange(TIME) %>%
  summarise(gc = list(SummarizeGrowth(TIME,ABS)), .groups = "drop") %>%
  mutate(
    vals = map(gc, ~ .x$vals)
  ) %>%
  mutate(
    param_df = map(vals, ~ as_tibble(.x[names(.x)]))
  ) %>%
  unnest_wider(param_df, names_repair = "unique")

Abs_export <- Abs_fits %>%
  select(where(~ !is.list(.)))
write.csv(Abs_export, "./data/12_invitro_growth/output/growthcurveR_data.csv", row.names = FALSE)


# Analysis mutants versus WT - Absorbance (ANOVA/Tukey, Welch/Games-Howell, Kruskal-Wallis/Dunn)
parameters <- c("k", "n0", "r", "t_mid", "t_gen", "auc_e", "auc_l")
pval_to_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "NS"
  )
}

Abs_fits_no_LB <- subset(Abs_fits,Abs_fits$STRAIN!="LB")

Abs_fits_results <- map_dfr(parameters, function(p) {
  formula <- as.formula(paste(p, "~ STRAIN"))
  model_data <- Abs_fits_no_LB %>% drop_na(!!sym(p))
  
  lm_model <- lm(formula, data = model_data)
  residuals <- resid(lm_model)
  shapiro_p <- tryCatch(shapiro.test(residuals)$p.value, error = function(e) NA)
  levene_p <- tryCatch(leveneTest(formula, data = model_data)$`Pr(>F)`[1], error = function(e) NA)
  anova_model <- tryCatch(aov(formula, data = model_data), error = function(e) NULL)
  anova_p <- tryCatch(summary(anova_model)[[1]][["Pr(>F)"]][1], error = function(e) NA)
  welch_p <- tryCatch(oneway.test(formula, data = model_data, var.equal = FALSE)$p.value, error = function(e) NA)
  kruskal_p <- tryCatch(kruskal.test(formula, data = model_data)$p.value, error = function(e) NA)
  
  if (!is.na(shapiro_p) && shapiro_p > 0.05) {
    if (!is.na(levene_p) && levene_p > 0.05 && !is.null(anova_model)) {
      tukey <- TukeyHSD(anova_model)
      results <- as_tibble(tukey[[1]], rownames = "contrast") %>%
        rename(p.adj = `p adj`, statistic = diff) %>%
        mutate(parameter = p,
               method = "Tukey",
               shapiro_p = shapiro_p,
               levene_p = levene_p,
               anova_p = anova_p,
               welch_p = welch_p,
               kruskal_p = kruskal_p,
               p.adj.signif = pval_to_stars(p.adj))
      
    } else {
      results <- model_data %>%
        games_howell_test(formula) %>%
        select(-.y.) %>%
        mutate(contrast = paste(group1, group2, sep = " - "),
               parameter = p,
               method = "Games-Howell",
               shapiro_p = shapiro_p,
               levene_p = levene_p,
               anova_p = anova_p,
               welch_p = welch_p,
               kruskal_p = kruskal_p,
               p.adj.signif = pval_to_stars(p.adj))
    }
  } else {
    dunn <- dunnTest(formula, data = model_data, method = "bonferroni")
    results <- dunn$res %>%
      rename(contrast = Comparison, statistic = Z, p.value = P.unadj, p.adj = P.adj) %>%
      mutate(parameter = p,
             method = "Dunn",
             shapiro_p = shapiro_p,
             levene_p = levene_p,
             anova_p = anova_p,
             welch_p = welch_p,
             kruskal_p = kruskal_p,
             p.adj.signif = pval_to_stars(p.adj))
  }
  
  results %>%
    select(contrast, parameter, shapiro_p, levene_p, method, 
           anova_p, welch_p, kruskal_p, p.adj, p.adj.signif, everything())
})

temp <- subset(Abs_fits_results, p.adj.signif != "NS")
write.csv(Abs_fits_results, "./data/12_invitro_growth/output/growthcurve_stats.csv", row.names = FALSE)