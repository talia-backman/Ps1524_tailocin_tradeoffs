# 06_plant_infections_luciferase.R
#
# Purpose:
#   Analyze bacterial fold-change (7 dpi) across strains/ecotypes and plot mean ± s.e.
#   Mixed-effects/saturated models per ecotype; pairwise vs WT with FDR correction.
#
# Input:
#   - data/06_plant_infections/luciferase_7dpi.csv
#
# Output:
#   - figures/06_plant_infections_luciferase.pdf
#   - Console output: model summaries, ANOVA, emmeans contrasts vs WT (FDR)
#
# Dependencies:
#   ggplot2, viridis, dplyr, lme4, lmerTest, emmeans
#
# Notes:
#   - Excludes buffer and OBC for main analyses (contamination/initial-lux issues).
#   - Uses lmer(fold_change_7dpi ~ strain + (1|batch)) if 'batch' exists, else lm().

library(ggplot2)
library(viridis)
library(lmerTest) # for p-values
library(emmeans)
library(dplyr)
library(lme4)

# read in data
dat <- read.csv("./data/06_plant_infections/luciferase_7dpi.csv")
# remove rows with NAs
dat <- dat[!(is.na(dat$ecotype) | dat$ecotype==""), ]
# remove buffer because it's inflated due to low initial lux values of plain buffer
dat <- subset(dat, strain != "buffer")
# remove OBC because we found contamination in one replicate and don't know how reliable the data is
dat <- subset(dat, strain != "OBC")

# get Col‑0 means for fold change
strain_ord <- dat %>% filter(ecotype == "Col0") %>%                  
  group_by(strain) %>% summarise(col0_mean = mean(fold_change_7dpi, na.rm = TRUE)) %>% 
  arrange(col0_mean) 
# make an ordered factor in dat
dat <- dat %>% mutate(strain_ord = factor(strain, levels = strain_ord$strain))
# Calculate the maximum y-value in the dataset for consistent placement of significance letters
y_max <- max(dat$fold_change_7dpi, na.rm = TRUE)
# Adjust the placement slightly above the max value
y_label_position <- y_max + 0.1
# plot, fold change lux
p <- ggplot(dat, aes(x = strain_ord, y = fold_change_7dpi, fill = strain_ord, colour = strain_ord)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", width = .75) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  geom_jitter(position = position_jitter(width = .05, height = 0)) +
  facet_wrap(~ ecotype) +
  scale_fill_viridis_d(option = "D") + scale_color_viridis_d(option = "D") +
  theme_bw(base_size = 20) + guides(fill = "none", colour = "none") +
  xlab("Strain") + ylab("Fold Change(Lux)") +
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = c("WT" = "WT",
                              "∆tagH" = expression(italic("∆tagH")),
                              "∆spsA" = expression(italic("∆spsA")),
                              "∆wfgD" = expression(italic("∆wfgD")),
                              "∆tagG" = expression(italic("∆tagG")),
                              "∆rmlC" = expression(italic("∆rmlC")),
                              "∆epsE" = expression(italic("∆epsE")),
                              "∆rmlC-FS" = expression(italic("∆rmlC-FS"))))
p

# cleaning for statistics per ecotype (vs WT)
# make sure strain is a factor
dat$strain <- factor(dat$strain)
# split by ecotype
col0  <- subset(dat, ecotype == "Col0")
eyach <- subset(dat, ecotype == "Eyach")
# set WT as reference
col0$strain  <- relevel(col0$strain,  ref = "WT")
eyach$strain <- relevel(eyach$strain, ref = "WT")

# do statistics
# ----- Col-0 -----
if ("batch" %in% names(col0)) {
  model_col0 <- lmer(fold_change_7dpi ~ strain + (1 | batch), data = col0)
} else {
  model_col0 <- lm(fold_change_7dpi ~ strain, data = col0)
}

print(summary(model_col0))
print(anova(model_col0))

# pairwise vs WT with FDR correction
emm_col0 <- emmeans(model_col0, "strain")
pairs_col0 <- contrast(emm_col0, method = "trt.vs.ctrl", ref = "WT", adjust = "fdr")
print(summary(pairs_col0, infer = TRUE))

# quick group means ± s.e.m.
col0_means <- col0 %>%
  group_by(strain) %>%
  summarise(
    n = sum(!is.na(fold_change_7dpi)),
    mean = mean(fold_change_7dpi, na.rm = TRUE),
    se = sd(fold_change_7dpi, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )
print(col0_means)

# ----- Eyach 1.5-2 -----
if ("batch" %in% names(eyach)) {
  model_eyach <- lmer(fold_change_7dpi ~ strain + (1 | batch), data = eyach)
} else {
  model_eyach <- lm(fold_change_7dpi ~ strain, data = eyach)
}
print(summary(model_eyach))
print(anova(model_eyach))

# pairwise vs WT with FDR correction
emm_eyach <- emmeans(model_eyach, "strain")
pairs_eyach <- contrast(emm_eyach, method = "trt.vs.ctrl", ref = "WT", adjust = "fdr")
print(summary(pairs_eyach, infer = TRUE))

# quick group means ± s.e.m.
eyach_means <- eyach %>%
  group_by(strain) %>%
  summarise(
    n = sum(!is.na(fold_change_7dpi)),
    mean = mean(fold_change_7dpi, na.rm = TRUE),
    se = sd(fold_change_7dpi, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )
print(eyach_means)

# re-plot with significance
dat <- dat %>%
  mutate(significance = case_when(
    strain == c("tagH","tagG","epsE","rmlC","wfgD","spsA","OBC") ~ "a",
    strain == c("WT") ~"b",
    TRUE ~ ""))
# plot, fold change lux
p <- ggplot(dat, aes(x = strain_ord, y = fold_change_7dpi, fill = strain_ord, colour = strain_ord)) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", width = .75) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  geom_jitter(position = position_jitter(width = .05, height = 0)) +
  facet_wrap(~ ecotype) +
  scale_fill_viridis_d(option = "D") + scale_color_viridis_d(option = "D") +
  theme_bw(base_size = 20) + guides(fill = "none", colour = "none") +
  xlab("Strain") + ylab("Fold Change(Lux)") +
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = c("WT" = "WT",
                              "tagH" = expression(italic("∆tagH")),
                              "spsA" = expression(italic("∆spsA")),
                              "wfgD" = expression(italic("∆wfgD")),
                              "tagG" = expression(italic("∆tagG")),
                              "rmlC" = expression(italic("∆rmlC")),
                              "epsE" = expression(italic("∆epsE")))) +
  geom_text(data = dat, aes(x = strain, y = y_label_position, label = significance),
            color = "black", size = 5) + theme(legend.position = "none")
p
# save plot
ggsave("./figures/06_plant_infections_luciferase.pdf", p, height = 5,width = 8)