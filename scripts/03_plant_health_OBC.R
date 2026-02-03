# 03_plant_health_OBC.R
#
# Purpose:
#   Compare plant health (mean green pixels) between OBC-present and OBC-absent groups.
#   Uses corrected OBC group assignments from OBC_presence.csv, removes strains without
#   explicit OBC calls, runs a Wilcoxon rank-sum test, and generates a boxplot figure.
#
# Input:
#   - data/input/plant_health_OBC_mean_pixels.csv
#   - data/input/plant_health_OBC_presence.csv
#
# Output:
#   - figures/plant_health_OBC.pdf
#   - data/output/plant_health_OBC_wilcoxon_summary.txt
#
# Dependencies:
#   dplyr, tidyr, ggplot2, viridis
#
# Notes:
#   - OBC_presence is treated as the source of truth (Absent/Present).
#   - Strains lacking an explicit OBC_presence assignment are removed.
#   - Non-parametric Wilcoxon test is used for robustness.

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# paths
in_mean_pixels <- "./data/input/plant_health_OBC_mean_pixels.csv"
in_obc_calls   <- "./data/input/plant_health_OBC_presence.csv"
out_fig   <- "./figures/plant_health_OBC.pdf"
out_dir   <- "./data/output/"
out_txt   <- file.path(out_dir, "plant_health_OBC_wilcoxon_summary.txt")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# load data
dat <- read.csv(in_mean_pixels, sep = " ")
hap <- read.csv(in_obc_calls)

hap <- hap %>%
  mutate(strain = gsub("plate", "p", strain)) %>%
  filter(!is.na(OBC_presence)) %>%
  mutate(OBC_presence = factor(OBC_presence, levels = c("Absent", "Present")))

obc_dat <- merge(dat, hap, by = "strain")

# filtering summary + stats
total_before <- nrow(obc_dat)

# if mean_pixel has NAs, drop them for test/plot
obc_filt <- obc_dat %>%
  filter(!is.na(OBC_presence)) %>%
  filter(!is.na(mean_pixel))

total_after <- nrow(obc_filt)

w_obc <- wilcox.test(mean_pixel ~ OBC_presence, data = obc_filt, exact = FALSE)

group_counts <- obc_filt %>%
  group_by(OBC_presence) %>%
  summarise(
    N = n(),
    strains = paste(strain, collapse = ", "),
    mean_mean_pixel = mean(mean_pixel, na.rm = TRUE),
    se_mean_pixel = sd(mean_pixel, na.rm = TRUE) / sqrt(sum(!is.na(mean_pixel))),
    .groups = "drop"
  )

# write summary txt
sink(out_txt)

cat("============================================================\n")
cat(" Plant health (mean green pixels) by OBC group\n")
cat("============================================================\n\n")

cat("Input files:\n")
cat(" - ", in_mean_pixels, "\n", sep = "")
cat(" - ", in_obc_calls, "\n\n", sep = "")

cat("Total rows after merge (before filtering): ", total_before, "\n", sep = "")
cat("Total rows after filtering (non-NA OBC + non-NA mean_pixel): ", total_after, "\n\n", sep = "")

cat("Group summary:\n")
print(group_counts)
cat("\n\n")

cat("Wilcoxon rank-sum test:\n")
print(w_obc)
cat("\n\n")

cat("Timestamp: ", as.character(Sys.time()), "\n", sep = "")

sink()

# plot: box plot (OBC− vs OBC+)
y_bar  <- max(obc_filt$mean_pixel, na.rm = TRUE) * 1.05
y_text <- y_bar * 1.01

p <- ggplot(obc_filt, aes(x = OBC_presence, y = mean_pixel, fill = OBC_presence)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = 16,
    outlier.size = 2,
    outlier.alpha = 0.9
  ) +
  scale_fill_viridis_d(option = "D") +
  theme_bw() +
  theme(
    text = element_text(size = 18),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15, lineheight = 0.9)
  ) +
  xlab("OBC Presence") +
  ylab("Mean green pixels") +
  geom_segment(aes(x = 1, xend = 2, y = y_bar, yend = y_bar), linewidth = 0.7) +
  annotate("text", x = 1.5, y = y_text, label = "**", size = 6, fontface = "bold")

ggsave(out_fig, p, width = 8, height = 6)
