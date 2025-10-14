# 04_tnseq_plant.R
#
# Purpose:
#   Process raw TnSeq read counts from *Pseudomonas syringae* p25C2 grown in planta
#   (Col-0 and Eyach 1.5-2 ecotypes) and merge resulting DESeq2-derived fitness
#   estimates with tailocin treatment data. These merged results form the basis
#   for identifying resistance trade-offs across host and tailocin conditions.
#
# Inputs:
#   - data/04_tnseq_scatter/input/plant_tnseq_counts.csv
#       Raw read counts per gene across host ecotypes and controls.
#   - data/04_tnseq_scatter/input/ortho.csv
#       Ortholog mapping between DC3000 and p25C2 gene IDs.
#   - data/04_tnseq_scatter/input/tailocin_tnseq.csv
#       TnSeq-derived logFC and p-values from tailocin exposure assays.
#
# Outputs:
#   - data/04_tnseq_scatter/output/merged_tnseq.csv
#       Unified table containing logFC and significance values for tailocin,
#       Col-0, and Eyach datasets, used for trade-off analyses and figures.
#
# Dependencies:
#   ggplot2, limma, DESeq2, tidyverse, dplyr, tibble
#
# Notes:
#   - Samples with <25,000 reads are filtered out.
#   - DC3000 IDs are cross-mapped to p25C2 using ortholog pairs.
#   - Only genes with non-missing tailocin data are retained for export.

library(ggplot2)
library(limma)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(tibble)

# Read in the counts table
all_exp <- read.table("./data/04_tnseq_scatter/input/plant_tnseq_counts.csv", header = TRUE, sep=",")

#there is inconsistency in naming. everything needs to be made lowercase
count_table_p <- t(all_exp[,12:dim(all_exp)[2]])
sample_order_p <- (all_exp[,1:11])
sample_order_p <- sample_order_p %>% mutate(plant = tolower(plant))

# Let's filter out all samples that have fewer than 25,000 reads
high_enough <- which(colSums(count_table_p, na.rm =TRUE)>=25000)
count_table_filter <-count_table_p[,high_enough]
sample_order_filter <- sample_order_p[high_enough,]

# Let's filter out all samples that are DC3000 treated
right_samples <- which(sample_order_filter$treatment!="dc3000" & sample_order_filter$experiment == "exp_0002" )
count_tab_filter <- count_table_filter[,right_samples]
count_tab_filter <-count_tab_filter[complete.cases(count_tab_filter), ]
samp_filter <- sample_order_filter[right_samples,]
samp_filter$plant <- factor(samp_filter$plant, levels = c("ctrl", "ey15_2", "col_0"))

#now with the filtered dataset, make deseq object
dds <- DESeqDataSetFromMatrix(countData = count_tab_filter, 
                              colData = samp_filter, design = ~ plant)
dds <- estimateSizeFactors(dds)

#now with the filtered dataset, make deseq object
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res_col <- as.data.frame(results(dds, contrast = c("plant","col_0", "ctrl" )))
res_col <- res_col[order(res_col$padj), ]
res_eyach <-  as.data.frame(results(dds, contrast = c("plant",  "ey15_2", "ctrl")))
colnames(res_eyach) <- paste0( colnames(res_eyach), "_eyach1025")
colnames(res_col) <- paste0( colnames(res_col), "_col1025")

# Give res_col an explicit ID column
res_col_tbl <- res_col %>% rownames_to_column("gene_id")
# Split that gene_id into p25 vs DC3000 helper columns
res_col_tbl <- res_col_tbl %>%
  mutate(p25_id    = if_else(grepl("^BJEIHDPM_", gene_id), gene_id, NA_character_),
         dc3000_id = if_else(grepl("^WP_",        gene_id), gene_id, NA_character_))

# read the DC3000 <-> p25.C2 mapping
# (file is tab-delimited and contains both WP_ and BJEIHDPM_ columns)
orth_map <- read_tsv("./data/04_tnseq_scatter/input/ortho.csv",show_col_types = FALSE) %>%
  transmute(dc3000_id = name_dc3000_new, p25_id = name_p25c2) %>% distinct()
# 2) Prepare the mapping and rename its p25 column to avoid suffix confusion
orth_map2 <- orth_map %>%
  rename(p25_id_map = p25_id)
# 3) Join and coalesce (no keep=TRUE; keeps a single dc3000_id)
res_col_filled <- res_col_tbl %>%
  left_join(orth_map2, by = "dc3000_id") %>%
  mutate(p25_id_final = dplyr::coalesce(p25_id, p25_id_map)) %>%
  select(-p25_id_map)
# sanity check
res_col_filled %>% summarise(n_total = n(), had_p25_already = sum(!is.na(p25_id)),          
                             filled_from_map = sum(is.na(p25_id) & !is.na(p25_id_final)), still_missing = sum(is.na(p25_id_final)))

# read in tailocin data
tailocin <- read.csv("./data/04_tnseq_scatter/input/tailocin_tnseq.csv", row.names = 1)
# make them all have an ID column to merge by
tailocin$ID <- sub(".*ID=([^;]+).*", "\\1", rownames(tailocin))
# ID tailocin columns
colnames(tailocin) <- paste0( colnames(tailocin), "_tailocin")
# merge tailocin with res_col_filled
merged_tradeoffs <- res_col_filled %>%
  left_join(tailocin, by = c("p25_id_final" = "ID_tailocin"))
# subset to only genes we have tailocin data for
merged_tradeoffs <- subset(merged_tradeoffs, logFC_tailocin != "NA")

# need to do the same for eyach dataframe
# Give res_eyach an explicit ID column
res_eyach_tbl <- res_eyach %>% rownames_to_column("gene_id") %>%
  mutate(p25_id    = if_else(grepl("^BJEIHDPM_", gene_id), gene_id, NA_character_), dc3000_id = if_else(grepl("^WP_", gene_id), gene_id, NA_character_))
# Fill using the ortholog map we already loaded
res_eyach_filled <- res_eyach_tbl %>%
  left_join(orth_map2, by = "dc3000_id") %>%
  mutate(p25_id_final = coalesce(p25_id, p25_id_map)) %>%
  select(-p25_id_map)

# merge eyach, col, and tailocin
# merge eyach with merged_tradeoffs
merged_tradeoffs_all <- merged_tradeoffs %>%
  left_join(res_eyach_filled, by = "p25_id_final")

# export data
write.csv(merged_tradeoffs_all,"./data/04_tnseq_scatter/output/merged_tnseq.csv")