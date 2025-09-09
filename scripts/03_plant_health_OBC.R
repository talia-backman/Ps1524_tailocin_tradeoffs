# 03_plant_health_OBC.R
#
# Purpose:
#   Compare plant health (mean green pixels) between OBC-present and OBC-absent groups.
#   Plot mean ± s.e. with jittered points; test OBC groups via Wilcoxon.
#
# Input:
#   - data/03_plant_health_OBC/mean_pixels.csv
#   - data/03_plant_health_OBC/OBC_presence.csv
#
# Output:
#   - figures/03_plant_health_OBC.pdf
#   - Wilcoxon test summary printed 
#
# Dependencies:
#   dplyr, ggplot2, viridis
#
# Notes:
#   - OBC factor levels standardized to c("Absent","Present").

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# read in data
dat <- read.csv("./data/03_plant_health_OBC/mean_pixels.csv", sep = " ")
hap <- read.csv("./data/03_plant_health_OBC/OBC_presence.csv")

# clean & merge
hap <- hap %>%
  mutate(strain = gsub("plate", "p", strain)) %>%
  # keep only strains where we have explicit OBC presence/absence
  filter(!is.na(OBC_presence)) %>%
  mutate(OBC_presence = factor(OBC_presence, levels = c("Absent", "Present")))

obc_dat <- merge(dat, hap, by = "strain")

# plot
# significance bar placement
y_bar <- max(obc_dat$mean_pixel, na.rm = TRUE) * 1.05
y_text <- y_bar * 1.01

set.seed(1)
p <- ggplot(obc_dat, aes(x = reorder(OBC_presence, mean_pixel, FUN = mean),
                         y = mean_pixel, fill = OBC_presence)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.75) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  geom_jitter(aes(colour = OBC_presence), width = 0.05, height = 0) +
  scale_fill_viridis_d(option = "D") +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "none") +
  xlab("OBC Presence") + ylab("Mean green pixels") +
  # significance bar (two groups: x=1 to x=2)
  geom_segment(aes(x = 1, xend = 2, y = y_bar, yend = y_bar)) +
  annotate("text", x = 1.5, y = y_text, label = expression(bold("**")), size = 6)

p
ggsave("./figures/03_plant_health_OBC.pdf", p, width = 8, height = 6)

# stats (Wilcoxon rank-sum, robust & consistent with prior script) 
w <- wilcox.test(mean_pixel ~ OBC_presence, data = obc_dat, exact = FALSE)

print(list(
  wilcox_test = w,
  group_sizes = table(obc_dat$OBC_presence),
  group_means = tapply(obc_dat$mean_pixel, obc_dat$OBC_presence, mean, na.rm = TRUE),
  group_se    = tapply(obc_dat$mean_pixel, obc_dat$OBC_presence,
                       function(x) sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x))))
))
