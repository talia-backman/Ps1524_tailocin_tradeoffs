# 10_aggregation.R
#
# Purpose:
#   Quantify strain-level aggregation phenotypes (R1–R3) and test whether O-antigen
#   mutants aggregate more frequently than WT.
#   - Plot stacked proportions of phenotypes per strain
#   - Global χ² test across strains
#   - Pairwise Fisher’s exact tests vs WT with FDR correction
#
# Inputs:
#   data/input/aggregation.csv
#     Required columns:
#       - strain
#       - phenotype  (e.g., "aggregation", "no aggregation", "clear")
#       - (optional) aggregation_binary (0/1 or text); if missing, derived from phenotype
#
# Outputs:
#   figures/aggregation.pdf
#   data/output/aggregation_chisq_overall.csv
#   data/output/aggregation_stats.csv
#
# Dependencies:
#   ggplot2, viridis, dplyr, tidyr, readr

library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(readr)

# ---------------- paths ----------------
in_file <- "data/input/aggregation.csv"
out_fig <- "figures/aggregation.pdf"
out_dir <- "data/output"
dir.create(dirname(out_fig), recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# read 
dat <- read.csv(in_file, check.names = FALSE)

# minimal clean: if aggregation_binary missing, derive from phenotype
if (!"aggregation_binary" %in% names(dat)) {
  dat$aggregation_binary <- ifelse(tolower(trimws(dat$phenotype)) == "aggregation", 1, 0)
}

# figure 
p <- ggplot(dat, aes(x = strain, fill = phenotype)) +
  geom_bar(position = "fill") +                                  
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(x = "Strain",
       y = "Proportion of wells",
       fill = "Phenotype",
       title = "Phenotype distribution by strain") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(out_fig, p, width = 7, height = 4, dpi = 300)

# stats
# Collapse to aggregation/non-aggregation (already binary)
tab <- with(dat, table(strain, aggregation_binary))

# Global chi-square across all strains
chi <- chisq.test(tab)
chisq_out <- data.frame(
  method    = chi$method,
  statistic = unname(chi$statistic),
  df        = unname(chi$parameter),
  p_value   = unname(chi$p.value),
  n_strains = nrow(tab),
  n_total   = nrow(dat)
)
write_csv(chisq_out, file.path(out_dir, "aggregation_chisq_overall.csv"))

# Pairwise Fisher vs WT (or p25C2 if that's your WT label)
baseline <- if ("WT" %in% rownames(tab)) "WT" else "p25C2"
if (!baseline %in% rownames(tab)) {
  baseline <- rownames(tab)[1]  # fallback so the script still runs
}

res <- lapply(setdiff(rownames(tab), baseline), function(s) {
  m <- rbind(tab[baseline, c("0","1")], tab[s, c("0","1")])
  dimnames(m) <- list(c(baseline, s), c("no_aggregation","aggregation"))
  ft <- fisher.test(m)
  data.frame(
    baseline   = baseline,
    strain     = s,
    odds_ratio = unname(ft$estimate),
    conf_low   = if (!is.null(ft$conf.int)) ft$conf.int[1] else NA_real_,
    conf_high  = if (!is.null(ft$conf.int)) ft$conf.int[2] else NA_real_,
    p_value    = ft$p.value
  )
})
res <- do.call(rbind, res)
res$p_adj <- p.adjust(res$p_value, method = "BH")
res <- res[order(res$p_adj, res$p_value), ]

# save
write_csv(res, file.path(out_dir, "aggregation_stats.csv"))