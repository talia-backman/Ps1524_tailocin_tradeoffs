# 02_HTF_OBC_sensitivity_summary.R
#
# Purpose:
#   Summarize the proportion of strains sensitive (S) or resistant (R) to tailocins
#   by presence/absence of the O-antigen biosynthesis cluster (OBC).
#   Converts long-format killing data into wide format, calculates phenotype counts,
#   and compares distributions between OBC-present and OBC-absent groups.
#
# Input:
#   - data/01_tailocin_killing_matrix/tailocin_killing_long.csv
#   - data/01_tailocin_killing_matrix/HTF_length_haplotype.csv
#
# Output:
#   - figures/02_HTF_OBC_sensitivity_summary.pdf
#   - Wilcoxon rank-sum test results printed to console
#
# Dependencies:
#   dplyr, tidyr, ggplot2, viridis
#
# Notes:
#   - OBC presence is derived from HTF length haplotypes
#   - Analyses exclude donor × tester combinations with fewer than 5 total counts
#   - Non-parametric Wilcoxon tests are used to compare OBC groups

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# read in data
dat <- read.csv("./data/01_tailocin_killing_matrix/tailocin_killing_long.csv")

# convert from long to wide form 
# add a replicate column by counting occurrences within each group
dat <- dat %>% 
  group_by(tester_strain, tailocin_donor_strain) %>%
  mutate(replicate = row_number()) %>% 
  ungroup()

# pivot wide with replicate numbers in column names
dat_wide <- dat %>% 
  pivot_wider(
    names_from = c(tailocin_donor_strain, replicate), 
    values_from = phenotype, 
    names_glue = "{tailocin_donor_strain}_tailocin_R{replicate}")

# calculate the counts of each phenotype per donor column
counts_list <- list()
for (col_name in colnames(dat_wide)[-1]) {  # exclude 'tester_strain' column
  counts <- table(dat_wide[[col_name]])
  count_df <- as.data.frame(counts)
  colnames(count_df) <- c("Phenotype", "Count")
  count_df$Column <- col_name
  count_df$Total <- sum(count_df$Count)
  counts_list[[col_name]] <- count_df
}
phenotype_counts_df <- bind_rows(counts_list)

# pivot to wide: one column for each phenotype (S = 1, R = 0)
phenotype_counts_wide <- phenotype_counts_df %>%
  pivot_wider(names_from = Phenotype, values_from = Count, 
    values_fill = list(Count = 0)) %>% select(Column, S, R, Total)

# add donor strain id by trimming suffix
phenotype_counts_wide$strain <- sub("_tailocin.*", "", phenotype_counts_wide$Column)

# read haplotype -> derive OBC presence/absence, then drop HTF entirely
hap <- read.csv("./data/01_tailocin_killing_matrix/HTF_length_haplotype.csv") %>%
  mutate(
    OBC = ifelse(HTF_length %in% c(1383, 1830), "OBC Absent", "OBC Present")
  ) %>% select(strain, OBC)  # keep only what's needed

# merge with counts
merged <- merge(phenotype_counts_wide, hap, by = "strain", all.x = TRUE)

# calculate percent S/R
merged <- merged %>%
  mutate(percent_S = S / Total,percent_R = R / Total)

# reshape for plotting and stats
merged_long <- merged %>%
  pivot_longer(
    cols = c(percent_S, percent_R),
    names_to = "Phenotype",
    values_to = "Percent"
  ) %>%
  filter(Total > 5) %>%              # keep inductions with enough counts
  filter(!is.na(OBC)) %>%            # drop rows without OBC annotation
  mutate(
    OBC = factor(OBC, levels = c("OBC Absent", "OBC Present"))
  )

# plot by OBC presence/absence
p <- ggplot(merged_long, aes(x = reorder(OBC, Percent, mean), y = Percent, fill = OBC)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.75) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  geom_jitter(aes(colour = OBC), width = 0.05, height = 0) +
  scale_fill_viridis_d(option = "D") +
  scale_color_viridis_d(option = "D") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "none") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap(~ Phenotype) +
  xlab("OBC Presence") +
  ylab("Percent of strains")

p
ggsave("./figures/02_HTF_OBC_sensitivity_summary.pdf", p, width = 8, height = 6)

# stats: compare OBC groups for each phenotype (non-parametric is fine/robust)
percent_R_data <- merged_long %>% filter(Phenotype == "percent_R")
percent_S_data <- merged_long %>% filter(Phenotype == "percent_S")

# Wilcoxon rank-sum tests (robust to non-normality & unequal n)
wR <- wilcox.test(Percent ~ OBC, data = percent_R_data, exact = FALSE)
wS <- wilcox.test(Percent ~ OBC, data = percent_S_data, exact = FALSE)

print(list(
  wilcox_percent_R = wR,
  wilcox_percent_S = wS,
  group_sizes_R = table(percent_R_data$OBC),
  group_sizes_S = table(percent_S_data$OBC)
))
