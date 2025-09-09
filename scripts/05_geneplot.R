# 05_geneplot.R
#
# Purpose:
#   Plot the focal OBC cluster gene map, colored by TnSeq results (Category).
#
# Input:
#   - data/05_geneplot/p25C2_tnseq_geneplot.csv
#
# Output:
#   - figures/05_geneplot.pdf
#
# Dependencies:
#   dplyr, ggplot2, viridis
#
# Notes:
#   - Ensure input has gene coordinates, strand, and an effect/score to map to color.

library(ggplot2)
library(gggenes) 
library(viridis)

# make color palette 
c25 <- c("black","#73D055FF")

# read gene plot data (a csv file with the following headers:
# molecule gene start end strand direction)
dat <- read.table("./data/05_geneplot/p25C2_tnseq_geneplot.csv", header = TRUE, sep=",")

# plot 
p <- ggplot(dat, aes(xmin = start, xmax = end, y = molecule, fill = Category)) +
  geom_gene_arrow() + labs(x = "Basepairs (bp)") + 
  ggtitle("O-Antigen biosynthesis gene cluster") +
  geom_text(aes(x = (start + end) / 2, y = as.numeric(as.factor(molecule)) + 0.08, 
                label = paste0("italic(", gene, ")")), 
            size = 3, angle = 30, hjust = 0, vjust = 0, parse = TRUE) + 
  scale_fill_manual(values = c25) + theme_minimal() +
  theme(axis.line = element_blank(),  # Removes axis lines
        axis.ticks = element_blank(), # Removes axis ticks
        axis.text.x = element_blank(),  # Removes axis text
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid = element_blank()) 
p

# save plot
ggsave("./figures/05_geneplot.pdf", p, width = 8, height = 6)
