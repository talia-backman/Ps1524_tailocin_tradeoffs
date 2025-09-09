# 04_tnseq_scatter.R
#
# Purpose:
#   Visualize TnSeq-derived plant fitness/phenotypes as a scatter (or similar summary).
#
# Input:
#   - data/04_tnseq_scatter/tailocin_plant_fitness.csv
#
# Output:
#   - figures/04_tnseq_scatter.pdf
#
# Dependencies:
#   dplyr, ggplot2, viridis
#
# Notes:
#   - Tailor axes/labels to fitness metric used in the CSV.

library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)

# read in TnSeq data
dat <- read.csv("./data/04_tnseq_scatter/tailocin_plant_fitness.csv", sep = ",")

# see how many genes there are per category
table(dat$Category)
# see how many tailocin only treatments were significant
table(dat$Significance, dat$Category)

# plot
p <- ggplot(dat, aes(x = log2FoldChange_Ey, y = logFC_tailocin)) +
  # Plot all other points
  geom_point(data = subset(dat, Category != "Resistance Trade-off"),
             aes(color = Significance_for_both, shape = Category, 
                 alpha = Significance_for_both)) +
  # Plot Resistance Trade-off separately with dark purple color
  geom_point(data = subset(dat, Category == "Resistance Trade-off"),
             aes(shape = Category, alpha = Significance_for_both), color = "#73D055FF") +
  scale_color_manual(values = c("Significant" = "black", "Not Significant" = "lightgrey")) +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.2), guide = "none") + 
  xlab("log Fold Change Eyach") + 
  ylab("log Fold Change Tailocin") +
  guides(color = guide_legend(title = "Significance")) + 
  theme_minimal()
p
# save 
ggsave("./figures/04_tnseq_scatter.pdf", p, width = 8, height = 6)

# add significant Resistance Trade-off gene labels
dat2_annotate <- subset(dat, Significance_for_both == "Significant" & Category == "Resistance Trade-off")
# Modify the plot to add labels
p2 <- ggplot(dat, aes(x = log2FoldChange_Ey, y = logFC_tailocin, 
                       color = Significance_for_both, shape = Category, 
                      alpha = Significance_for_both)) +
  geom_point() +  
  scale_color_manual(values = c("Significant" = "black", "Not Significant" = "lightgrey")) +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.4), guide = "none") + 
  xlab("log Fold Change Eyach") + 
  ylab("log Fold Change Tailocin") +
  guides(color = guide_legend(title = "Significance")) + 
  theme_minimal() +
  geom_text_repel(data = dat2_annotate, aes(label = name_p25c2_new), size = 2, box.padding = 0.1)
p2 # for annotating genes in paper
