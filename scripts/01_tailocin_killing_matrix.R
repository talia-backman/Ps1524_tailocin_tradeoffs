# 01_tailocin_killing_matrix.R
# 
# Purpose:
#   Visualize the tailocin killing matrix alongside a phylogenetic tree of tester strains.
#   Adds tip labels colored by HTF length haplotype and a heatmap showing killing outcomes.
#
# Input:
#   - data/01_tailocin_killing_matrix/tailocin_killing_matrix.csv
#   - data/01_tailocin_killing_matrix/HTF_length_haplotype.csv
#   - data/01_tailocin_killing_matrix/ps_1524_uncollapsed_5_2018.nwk
#
# Output:
#   - figures/01_tailocin_killing_matrix.pdf
#
# Dependencies:
#   ggplot2, ape, ggtree, dplyr, viridis, RColorBrewer
#
# Notes:
#   - Subsets the tree to include only tester strains present in the killing matrix.
#   - Columns are reordered to match haplotype groups before plotting the heatmap.

library(ggplot2)
library(ape) 
library(ggtree) 
library(dplyr) 
library(viridis) 
library(RColorBrewer) 

# read in data and prep it for plotting
dat <- read.csv("./data/01_tailocin_killing_matrix/tailocin_killing_matrix.csv")
htf_dat <- read.csv("./data/01_tailocin_killing_matrix/HTF_length_haplotype.csv")
# merge data
dat <- merge(dat, htf_dat, by.x = "X", by.y = "strain")
# read in tree
tree <- read.tree("./data/01_tailocin_killing_matrix/ps_1524_uncollapsed_5_2018.nwk")
# subset tree to only the tester strains used in dat
uneeded <- which(!(tree$tip.label %in% dat$X))
tree <- drop.tip(tree, uneeded)
length(tree$tip.label) # sanity check
length(unique(dat$X)) # sanity check
# color phylogeny tip labels
tree_data <- dat %>% select(X, HTF_length) %>%
  mutate(HTF_length = as.factor(HTF_length))
# Define custom colors for haplotypes
tree_palette <- c(
  "1830" = "#e8b9d8",  # light pink
  "1383" = "#4a76b7",  # blue
  "1803" = "#8c0027",  # dark red
  "1245" = "#e6b625")   # mustard yellow
# make a copy of the data for heatmap, keep dat (with HTF_length) for tip colors
dat_heat <- dat
rownames(dat_heat) <- dat_heat[,1]
dat_heat <- dat_heat[ ,c(2:8)]  # 7 assay columns after the ID column
# convert to matrix for plotting
dat_matrix <- as.matrix(dat_heat)
# define the correct order for plotting
correct_order <- c("p7.G11","p23.B2", # 1830bp
                   "p25.A12","p1.G2", #1383bp
                   "p5.C3","p21.F9", #1245bp
                   "p23.B8") # 1803bp
# reorder columns
dat_matrix <- dat_matrix[, correct_order]

# plot
# Create the base tree plot with colored tip labels
p <- ggtree(tree) %<+% dat + 
  geom_tiplab(aes(color = as.factor(HTF_length)), size = 2) +  
  scale_color_manual(values = tree_palette, name = "HTF Length (Tester)")
p
# Add heatmap with colored column annotations
p3 <- gheatmap(p, dat_matrix, colnames = TRUE, legend_title = "Results", 
               offset = .05, color = "black", colnames_position = "bottom", 
               font.size = 3, width = 5) + scale_fill_viridis(discrete = FALSE) 
p3

# save figure as pdf
ggsave("./figures/01_tailocin_killing_matrix.pdf", p3, width = 12, height = 10)
