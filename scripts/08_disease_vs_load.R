# 08_disease_vs_load.R
#
# Purpose:
#   Relate plant health (proportion Healthy) to bacterial load (mean fold-change, 7 dpi).
#   Builds %Healthy by strain/ecotype from 07, joins to 06 means, fits per-ecotype LMs,
#   and generates an annotated correlation plot.
#
# Input:
#   - data/07_plant_disease/plant_disease.csv
#   - data/07_plant_disease/contam.csv
#   - data/07_plant_disease/combined_plate_layouts.csv
#   - data/06_plant_infections/luciferase_7dpi.csv
#
# Output:
#   - figures/08_disease_vs_load.pdf
#   - results/08_disease_vs_load/per_ecotype_linear_models.csv
#
# Dependencies:
#   dplyr, tidyr, ggplot2, viridis, broom
#
# Notes:
#   - Excludes buffer and OBC to match main analyses.
#   - Reports per-ecotype R^2 and p-values; pooled model printed for reference.

library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(broom)

# Build %Healthy per strain/ecotype from 07 inputs
# read in qualitative plant data
dat <- read.csv("./data/07_plant_disease/plant_disease.csv")

# collapse 4 categories -> 2
dat$Phenotype[dat$Phenotype == "B"] <- "Healthy"
dat$Phenotype[dat$Phenotype == "S"] <- "Healthy"
dat$Phenotype[dat$Phenotype == "X"] <- "Diseased/Dead"
dat$Phenotype[dat$Phenotype == "D"] <- "Diseased/Dead"

# read in contamination rows to remove
contam <- read.csv("./data/07_plant_disease/contam.csv")

# drop contaminated wells
dat_filtered <- dat %>%
  anti_join(contam, by = c("Well", "File_Name"))

# plate layout (to get strain/ecotype keys)
layout_plate <- read.csv("./data/07_plant_disease/combined_plate_layouts.csv", header = TRUE)

# merge by Well + file name (file_name in layout vs File_Name in scoring)
merged_data <- merge(layout_plate, dat_filtered,
                     by.x = c("Well", "file_name"),
                     by.y = c("Well", "File_Name"),
                     all.x = TRUE)

# remove NAs and rmlC-FS 
merged_data <- na.omit(merged_data)
merged_data <- subset(merged_data, strain != "RmlC-FS")

# derive ecotype from file_name
merged_data$ecotype <- sub("^(Col0|Eyach)_.*", "\\1", merged_data$file_name)

# % Healthy per strain/ecotype
plot_data <- merged_data %>%
  group_by(strain, ecotype) %>%
  summarise(
    total   = n(),
    X_count = sum(Phenotype == "Healthy"),
    .groups = "drop"
  ) %>%
  mutate(prop_X = X_count / total)

# load in 06 infection-load means and join
lux_dat <- read.csv("./data/06_plant_infections/luciferase_7dpi.csv")

mean_lux <- lux_dat %>%
  filter(!is.na(ecotype), ecotype != "") %>%
  group_by(strain, ecotype) %>%
  summarise(
    mean_fold_change = mean(fold_change_7dpi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  na.omit()

# join and drop buffer/OBC to match earlier analyses
combined_dat <- inner_join(mean_lux, plot_data, by = c("strain", "ecotype")) %>%
  filter(!strain %in% c("buffer", "OBC"))

# stats: per-ecotype linear models (prop_X ~ mean_fold_change)
# per-ecotype models
model_stats <- combined_dat %>%
  group_by(ecotype) %>%
  do({
    mod <- lm(prop_X ~ mean_fold_change, data = .)
    tibble(
      r_squared = summary(mod)$r.squared,
      slope     = coef(mod)[["mean_fold_change"]],
      p_value   = summary(mod)$coefficients["mean_fold_change", "Pr(>|t|)"]
    )
  }) %>% ungroup()
print(model_stats)

# pooled model (optional)
model_pooled <- lm(prop_X ~ mean_fold_change + ecotype, data = combined_dat)
print(summary(model_pooled))

# plot: proportion Healthy vs mean fold-change, per ecotype trend lines
# annotation strings
lab_col0  <- model_stats %>% filter(ecotype == "Col0")  %>%
  transmute(lbl = paste0("Col-0: R^2 = ", round(r_squared, 2),
                         ", p = ", signif(p_value, 2))) %>% pull(lbl)
lab_eyach <- model_stats %>% filter(ecotype == "Eyach") %>%
  transmute(lbl = paste0("Eyach: R^2 = ", round(r_squared, 2),
                         ", p = ", signif(p_value, 2))) %>% pull(lbl)

# make the plot
p <- ggplot(combined_dat, aes(x = mean_fold_change, y = prop_X, color = ecotype, label = strain)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, aes(group = ecotype), linetype = "dashed", size = 1) +
  geom_text(vjust = -0.6, size = 3) +
  scale_color_viridis_d(option = "D") +
  theme_bw(base_size = 20) +
  xlab("Mean bacterial load (fold change, 7 dpi)") +
  ylab("Proportion healthy") +
  ggtitle("Plant health vs. bacterial load across ecotypes") +
  theme(legend.position = "top")

# place annotations in the upper right of the current data range
x_pos <- max(combined_dat$mean_fold_change, na.rm = TRUE)
y_pos1 <- max(combined_dat$prop_X, na.rm = TRUE)
y_pos2 <- y_pos1 - 0.08

p <- p +
  annotate("text", x = x_pos, y = y_pos1, hjust = 1, label = lab_col0,  size = 5, color = "black") +
  annotate("text", x = x_pos, y = y_pos2, hjust = 1, label = lab_eyach, size = 5, color = "black")

p
# save plot
ggsave("./figures/08_disease_vs_load.pdf", p, width = 8, height = 6)