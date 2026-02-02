# 02_HTF_OBC_sensitivity_summary.R
#
# Purpose:
#   Summarize tailocin killing outcomes as the % of tester strains killed (S) vs not killed (R)
#   for each tailocin donor strain replicate, and compare distributions between OBC+ vs OBC-
#   donor groups. Produces boxplots and Wilcoxon rank-sum tests.
#
# Inputs:
#   - data/01_tailocin_killing_matrix/tailocin_killing_long.csv
#       Long-format killing assay results with (at minimum) columns:
#         tester_strain, tailocin_donor_strain, phenotype
#       where phenotype is "S" (killed) or "R" (not killed).
#   - data/01_tailocin_killing_matrix/HTF_length_haplotype.csv
#       Mapping of donor strain -> HTF_length with (at minimum) columns:
#         strain, HTF_length
#
# Outputs:
#   - figures/02_HTF_OBC_sensitivity_summary.pdf
#       Two-page PDF:
#         Page 1: % of strains killed (percent_S) by OBC group
#         Page 2: % of strains not killed (percent_R) by OBC group
#   - figures/02_HTF_OBC_sensitivity_summary_stats.txt
#       Text summary including group sizes and Wilcoxon test results.
#
# Dependencies:
#   dplyr, tidyr, ggplot2, viridis
#
# Notes:
#   - OBC group is derived from HTF_length haplotypes:
#       HTF_length in {1383, 1830} => "OBC-"  (OBC Absent)
#       otherwise                  => "OBC+"  (OBC Present)
#   - Each donor strain replicate is treated as one observation:
#       percent_S / percent_R computed across all tester strains for that donor replicate.
#   - Inductions/replicates with Total <= 5 are excluded (insufficient tester counts).
#   - Wilcoxon rank-sum tests compare OBC- vs OBC+ distributions for percent_S and percent_R.

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# paths
in_killing <- "./data/01_tailocin_killing_matrix/tailocin_killing_long.csv"
in_hap     <- "./data/01_tailocin_killing_matrix/HTF_length_haplotype.csv"
out_pdf    <- "./figures/02_HTF_OBC_sensitivity_summary.pdf"
out_txt    <- "./data/02_HTF_OBC_sensitivity_summary/02_HTF_OBC_sensitivity_summary_stats.txt"

# load + format data
dat <- read.csv(in_killing)

# add replicate index within each tester × donor group
dat <- dat %>%
  group_by(tester_strain, tailocin_donor_strain) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

# pivot to wide so each donor replicate becomes its own column of phenotypes across testers
dat_wide <- dat %>%
  pivot_wider(
    names_from = c(tailocin_donor_strain, replicate),
    values_from = phenotype,
    names_glue = "{tailocin_donor_strain}_R{replicate}"
  )

# count S/R per donor replicate column
count_list <- list()
for (col_name in colnames(dat_wide)[-1]) {  # exclude tester_strain
  tbl <- table(dat_wide[[col_name]])
  df <- as.data.frame(tbl)
  colnames(df) <- c("Phenotype", "Count")
  df$Column <- col_name
  df$Total <- sum(df$Count)
  count_list[[col_name]] <- df
}
counts_df <- bind_rows(count_list)

counts_wide <- counts_df %>%
  pivot_wider(names_from = Phenotype, values_from = Count, values_fill = 0) %>%
  select(Column, S, R, Total)

# extract donor strain name from column label
counts_wide$strain <- sub("_R.*", "", counts_wide$Column)

# load haplotype mapping and derive OBC group
hap <- read.csv(in_hap) %>%
  mutate(OBC = ifelse(HTF_length %in% c(1383, 1830), "OBC-", "OBC+")) %>%
  select(strain, OBC)

# merge + compute percentages; filter to annotated donors and sufficient totals
merged <- merge(counts_wide, hap, by = "strain", all.x = TRUE) %>%
  filter(!is.na(OBC), Total > 5) %>%
  mutate(
    percent_S = S / Total * 100,
    percent_R = R / Total * 100,
    OBC = factor(OBC, levels = c("OBC-", "OBC+"))
  )

# two-line x labels for plotting
merged$OBC_label <- factor(
  merged$OBC,
  levels = c("OBC-", "OBC+"),
  labels = c("Tailocins from\nOBC- strains", "Tailocins from\nOBC+ strains")
)

# STATS (WILCOXON)
percent_S_data <- merged %>% select(strain, Column, OBC, OBC_label, percent_S)
percent_R_data <- merged %>% select(strain, Column, OBC, OBC_label, percent_R)

wS <- wilcox.test(percent_S ~ OBC, data = percent_S_data, exact = FALSE)
wR <- wilcox.test(percent_R ~ OBC, data = percent_R_data, exact = FALSE)

group_counts <- merged %>%
  group_by(OBC) %>%
  summarise(
    N = n(),
    donor_replicates = paste(Column, collapse = ", "),
    .groups = "drop"
  )

# helper for significance annotation
p_to_stars <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("NS")
}

stars_S <- p_to_stars(wS$p.value)
stars_R <- p_to_stars(wR$p.value)

# Write summary text file
sink(out_txt)
cat("============================================\n")
cat(" OBC Tailocin Sensitivity Summary\n")
cat("============================================\n\n")

cat("Inputs:\n")
cat(" - ", in_killing, "\n", sep = "")
cat(" - ", in_hap, "\n\n", sep = "")

cat("Filtering:\n")
cat(" - Kept donor replicates with Total > 5 tester outcomes\n")
cat(" - Dropped rows lacking OBC annotation\n\n")

cat("Group sample counts (each row = one donor replicate column):\n")
print(group_counts)
cat("\n\n")

cat("Wilcoxon test — Percent Sensitive (percent_S):\n")
print(wS)
cat("\nSignificance: ", stars_S, "\n\n", sep = "")

cat("Wilcoxon test — Percent Resistant (percent_R):\n")
print(wR)
cat("\nSignificance: ", stars_R, "\n\n", sep = "")

cat("Timestamp: ", as.character(Sys.time()), "\n", sep = "")
sink()

# PLOTTING
make_plot <- function(df, yvar, ylab_txt, stars_label) {
  y_max <- max(df[[yvar]], na.rm = TRUE)
  ybar <- y_max * 1.05
  ytxt <- y_max * 1.10
  
  ggplot(df, aes(x = OBC_label, y = .data[[yvar]], fill = OBC_label)) +
    geom_boxplot(
      alpha = 0.7,
      outlier.shape = 16,
      outlier.size = 2,
      outlier.alpha = 0.9
    ) +
    geom_jitter(width = 0.10, alpha = 0.8, size = 2) +
    scale_fill_viridis_d(option = "D") +
    theme_bw() +
    theme(
      text = element_text(size = 18),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 15, lineheight = 0.9)
    ) +
    ylab(ylab_txt) +
    coord_cartesian(ylim = c(0, ytxt * 1.05)) +
    geom_segment(aes(x = 1, xend = 2, y = ybar, yend = ybar), linewidth = 0.8) +
    annotate("text", x = 1.5, y = ytxt, label = stars_label, size = 6, fontface = "bold")
}

pS <- make_plot(merged, "percent_S", "% of strains killed", stars_S)
pR <- make_plot(merged, "percent_R", "% of strains not killed", stars_R)

# save as PDF
pdf(out_pdf, width = 8, height = 6)
print(pS)
print(pR)
dev.off()
