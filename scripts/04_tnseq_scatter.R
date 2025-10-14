# 04_tnseq_scatter.R
#
# Purpose:
#   Visualize TnSeq-derived phenotypes and highlight "resistance trade-off" genes
#   (tailocin ↑, plant fitness ↓) across Arabidopsis ecotypes.
#
# Input:
#   - data/04_tnseq_scatter/input/merged_tnseq.csv
#
# Outputs:
#   - figures/04_tnseq_scatter.pdf
#   - data/04_tnseq_scatter/output/manuscript_table.csv
#
# Dependencies:
#   ggplot2, dplyr, viridis, ggrepel, readr

library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(stringr)
library(ggnewscale)
library(readr)

# read in TnSeq data
dat <- read.csv("./data/04_tnseq_scatter/input/merged_tnseq.csv")

# add significance column for plotting
dat <- dat %>% mutate(Significance = ifelse(pvalue_eyach1025 < 0.05 & P.Value_tailocin < 0.05 & pvalue_col1025 < 0.05, 
                                            "Significant", "Not Significant"))

# add category column for plotting
dat <- dat %>%
  mutate(Category = case_when(
    log2FoldChange_eyach1025 < 0 & log2FoldChange_col1025 < 0 & logFC_tailocin > 0 ~ "Resistance Trade-off both ecotypes",
    log2FoldChange_eyach1025 > 0 & log2FoldChange_col1025 > 0 & logFC_tailocin > 0 ~ "Mutual Benefit",
    log2FoldChange_eyach1025 < 0 & log2FoldChange_col1025 < 0 & logFC_tailocin < 0 ~ "Mutual Detriment",
    log2FoldChange_eyach1025 > 0 & log2FoldChange_col1025 > 0 & logFC_tailocin < 0 ~ "Non-Focal Trade-off"))
# see how many genes there are per category
table(dat$Category)
# see how many tailocin only treatments were significant
table(dat$Significance, dat$Category) # 14 but duplicate row..
# subset resistance trade off (RT) genes
RT <- subset(dat, Category == "Resistance Trade-off both ecotypes") 
RT_sig <- subset(RT, Significance == "Significant")

# define RT_set by effect directions (tailocin > 0 & plant < 0)
dat2 <- dat %>%
  mutate(
    RT_flag_col  = logFC_tailocin > 0 & log2FoldChange_col1025  < 0,
    RT_flag_ey   = logFC_tailocin > 0 & log2FoldChange_eyach1025 < 0,
    RT_set = case_when(RT_flag_col & RT_flag_ey ~ "Both",
                       RT_flag_col & !RT_flag_ey ~ "Col-0 only",!RT_flag_col & RT_flag_ey ~ "Eyach only",TRUE ~ "None"))
# define RT sets (tailocin > 0 & plant < 0 with p<0.05)
alpha_thr <- 0.05
dat2 <- dat2 %>% mutate(Significance_RT = case_when(
  RT_set == "Both"       ~ ifelse(P.Value_tailocin < alpha_thr &
                                    pvalue_col1025  < alpha_thr &
                                    pvalue_eyach1025 < alpha_thr, "Significant", "Not Significant"),
  RT_set == "Col-0 only" ~ ifelse(P.Value_tailocin < alpha_thr &
                                    pvalue_col1025  < alpha_thr, "Significant", "Not Significant"),
  RT_set == "Eyach only" ~ ifelse(P.Value_tailocin < alpha_thr &
                                    pvalue_eyach1025 < alpha_thr, "Significant", "Not Significant"),
  TRUE ~ NA_character_))

# eyach scatterplot with colored crosses for RT sets 
rt_cols <- c("Both" = "#73D055FF","Col-0 only" = "#440154FF","Eyach only" = "#31688EFF")

# Only annotate selected genes
dat2_annotate <- subset(dat2, gene_or_product %in% c("wfgD", "rmlC_1", "tagG_2", "tagH_2", "epsE_4"))
# for some reason rmlC_1 duplicates, so remove duplicate row
dat2_annotate <- dat2_annotate %>% distinct(gene_or_product, .keep_all = TRUE)
# plot!
p3 <- ggplot(dat2, aes(x = log2FoldChange_eyach1025, y = logFC_tailocin)) +
  geom_point(data = subset(dat2, Category != "Resistance Trade-off"),
             aes(color = Significance, shape = Category, alpha = Significance), size = 1.4) +
  geom_point(data = subset(dat2, RT_set != "None"),
             aes(color = RT_set, shape = Category, alpha = Significance_RT), size = 1.4, stroke = 1.2) +
  scale_color_manual(values = c("Significant" = "black", "Not Significant" = "lightgrey", rt_cols),
                     breaks = c("Not Significant", "Significant"), name = "Significance") +
  scale_alpha_manual(values = c("Significant" = 1, "Not Significant" = 0.2), name = "Significance") +
  xlab("log Fold Change Eyach") +
  ylab("log Fold Change Tailocin") +
  theme_minimal() +
  geom_text_repel(data = dat2_annotate, aes(label = gene_or_product),
                  size = 2, box.padding = 0.1)
p3
# save
ggsave("./figures/04_tnseq_scatter.pdf", p3, width = 8, height = 4, dpi = 300)


#################### for generating the tnseq significant resistance trade off genes table #################### 
library(readr)
library(dplyr)

# load in data
df <- read_csv("./data/04_tnseq_scatter/input/merged_tnseq.csv", show_col_types = FALSE)
# subset to only genes in resistance trade off category
df <- subset(df, logFC_tailocin > 0)
df <- subset(df, P.Value_tailocin < 0.05) # 71 obs. 
# now make significance column
df <- df %>%
  mutate(significance = case_when(
    (pvalue_eyach1025 < 0.05 & pvalue_col1025 < 0.05) ~ "both",
    (pvalue_col1025 < 0.05 & pvalue_eyach1025 > 0.05) ~ "col0 only",
    (pvalue_eyach1025 < 0.05 & pvalue_col1025 > 0.05) ~ "eyach only",
    TRUE ~ "not significant"))
# subset to only significant interactions
df <- subset(df, significance != "not significant") # 36 obs.

# add logFC column for resistance trade off genes
df <- df %>% mutate(logFC_category = case_when(
  (log2FoldChange_eyach1025 < 0 & log2FoldChange_col1025 < 0) ~ "logFC_both",
  (log2FoldChange_col1025 < 0 & log2FoldChange_eyach1025 > 0) ~ "logFC_col0 only",
  (log2FoldChange_eyach1025 < 0 & log2FoldChange_col1025 > 0) ~ "logFC_eyach only",
  TRUE ~ NA_character_))
# subset rows with NA in logFC column out of df
df <- subset(df, logFC_category != "NA")

# see if there are duplicate rows..
df <- df[!duplicated(df), ] # none

# counts
table(df$logFC_category) # 24
# see how many tailocin only treatments were significant
table(df$significance) # 24 

# write csv
write_csv(df, "./data/04_tnseq_scatter/output/manuscript_table.csv")