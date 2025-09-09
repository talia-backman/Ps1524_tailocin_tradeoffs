# 07_plant_disease.R
#
# Purpose:
#   Produce stacked bars of disease outcomes (Healthy vs Diseased/Dead) per strain/ecotype,
#   and test whether mutants differ from WT in Healthy proportion (logistic GLM + FDR).
#
# Input:
#   - data/07_plant_disease/plant_disease.csv
#   - data/07_plant_disease/contam.csv
#   - data/07_plant_disease/combined_plate_layouts.csv
#
# Output:
#   - figures/07_plant_disease.pdf
#   - Console output: Dunnett-style comparisons vs WT (GLM, FDR-adjusted)
#
# Dependencies:
#   dplyr, tidyr, ggplot2, viridis, multcomp
#
# Notes:
#   - Collapses four phenotype codes (B,S,X,D) into two: Healthy vs Diseased/Dead.
#   - Removes contaminated wells via anti_join(Well, File_Name).
#   - Excludes OBC for main plots/tests (data quality considerations).

library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr) 
library(multcomp)

# read in qualitative plant data
# I manually re-named the split plate (Col0-Eyach_R3) to match contam data frame
dat <- read.csv("./data/07_plant_disease/plant_disease.csv")
# replace 4 groups with only 2
dat$Phenotype[dat$Phenotype == "B"] <- "Healthy"
dat$Phenotype[dat$Phenotype == "S"] <- "Healthy"
dat$Phenotype[dat$Phenotype == "X"] <- "Diseased/Dead"
dat$Phenotype[dat$Phenotype == "D"] <- "Diseased/Dead"

# read in contamination data (rows to remove from dat)
contam <- read.csv("./data/07_plant_disease/contam.csv")

# Remove rows from dat where both 'Well' and 'file_name' match those in contam
dat_filtered <- dat %>%
  anti_join(contam, by = c("Well", "File_Name"))
# see if there are duplicate rows in contam
contam %>% group_by(Well, File_Name) %>% filter(n() > 1) # no duplicates
# read in plate layout
layout_plate <- read.csv("./data/07_plant_disease/combined_plate_layouts.csv", header = T)
# Merge the two dataframes based on exact matches of Well and file_name <-> Batch
merged_data <- merge(layout_plate, dat_filtered, 
                     by.x = c("Well", "file_name"), 
                     by.y = c("Well", "File_Name"), all.x = TRUE)
# remove NAs
merged_data <- na.omit(merged_data)
# remove rmlC-FS
merged_data <- subset(merged_data, strain != "RmlC-FS")
# Create the new column 'ecotype' based on the first part of 'file_name'
merged_data$ecotype <- sub("^(Col0|Eyach)_.*", "\\1", merged_data$file_name)

# First, count how many Phenotype X per strain and total counts
plot_data <- merged_data %>%
  group_by(strain, ecotype) %>%
  summarise(
    total = n(),
    X_count = sum(Phenotype == "Healthy")
  ) %>%
  mutate(prop_X = X_count / total)

# Now join back with full count data (for stacked bar plot)
full_plot_data <- merged_data %>%
  count(strain, Phenotype, ecotype) %>%
  left_join(plot_data %>% select(strain, ecotype, prop_X), by = c("strain", "ecotype")) %>%
  mutate(strain = reorder(strain, prop_X))  # <-- reorder strains based on proportion of X

# remove OBC from data
full_plot_data <- subset(full_plot_data, strain != "OBC")
# Now plot
p <- ggplot(full_plot_data, aes(x = strain, y = n, fill = Phenotype)) +
  geom_bar(stat = "identity",position = "fill") +
  facet_wrap(~ ecotype) +
  theme_bw() + theme(text = element_text(size = 20)) +
  xlab("Strain") +
  ylab("Percent") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis_d(option = "D") +
  scale_x_discrete(guide = guide_axis(angle = 90),
                   labels = c("WT" = "WT",
                              "tagH" = expression(italic("∆tagH")),
                              "spsA" = expression(italic("∆spsA")),
                              "wfgD" = expression(italic("∆wfgD")),
                              "tagG" = expression(italic("∆tagG")),
                              "rmlC" = expression(italic("∆rmlC")),
                              "epsE" = expression(italic("∆epsE"))))
p
# save as pdf
ggsave("./figures/07_plant_disease.pdf", p, width = 7, height = 4)

# stats
# Make binary outcome: Healthy = 1, Diseased/Dead = 0
binom_dat <- full_plot_data %>%
  mutate(outcome = ifelse(Phenotype == "Healthy", 1, 0))

# Expand counts (n = number of plants) into individual rows
glm_dat <- binom_dat %>%
  uncount(weights = n)  # repeats each row 'n' times

# Fit logistic regression separately per ecotype
col0  <- subset(glm_dat, ecotype == "Col0")
eyach <- subset(glm_dat, ecotype == "Eyach")

# WT as reference
col0$strain  <- relevel(factor(col0$strain), ref = "WT")
eyach$strain <- relevel(factor(eyach$strain), ref = "WT")

# Logistic regressions
model_col0  <- glm(outcome ~ strain, family = binomial, data = col0)
model_eyach <- glm(outcome ~ strain, family = binomial, data = eyach)

summary(model_col0)
summary(model_eyach)

# Multiple comparisons vs WT (Dunnett-style)
comp_col0 <- glht(model_col0, linfct = mcp(strain = "Dunnett"))
comp_eyach <- glht(model_eyach, linfct = mcp(strain = "Dunnett"))

summary(comp_col0, test = adjusted("fdr"))
summary(comp_eyach, test = adjusted("fdr"))

summary(comp_col0, test = adjusted("fdr"))
summary(comp_eyach, test = adjusted("fdr"))