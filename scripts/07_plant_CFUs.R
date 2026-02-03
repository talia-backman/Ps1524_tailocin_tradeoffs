# 13_plant_CFUs.R
#
# Purpose:
#   Analyze in planta CFU data from synthesis infections (WT vs wfgD across inoculum doses)
#   in two Arabidopsis ecotypes (Eyach 1.5-2 and Col-0) at 2 and 5 dpi. Generates per-ecotype,
#   per-timepoint boxplots and runs pairwise WT vs wfgD comparisons within each inoculum dose.
#
# Input:
#   - data/input/plant_CFUs_Synthesis.csv
#
# Outputs:
#   - figures/13_Eyach_2dpi.pdf
#   - figures/13_Eyach_5dpi.pdf
#   - figures/13_Col0_2dpi.pdf
#   - figures/13_Col0_5dpi.pdf
#   - data/output/plant_CFUs_stats.csv
#
# Notes:
#   - MODA is constructed as TREATMENT_CONCENTRATION (e.g., WT_10^6).
#   - Plots use LOG10.CFU.LEAF. as the response variable and fix y-limits to 0–8.
#   - Statistical testing is adaptive per concentration (t-test/Welch/Wilcoxon) based on normality/variance checks.
#
# Dependencies:
#   tidyr, tidyverse, dplyr, ggplot2, purrr, broom, FSA, carData, car, rstatix, dunn.test, stringr
#

#Import useful libraries
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(purrr)
library(broom)
library(FSA)
library(carData)
library(car)
library(rstatix)
library(dunn.test)
library(stringr)

# Import transform and organize data
data <- read.csv("./data/input/plant_CFUs_Synthesis.csv")
data$DPI <- as.numeric(data$DPI)
data$CFU.LEAF <- as.numeric(data$CFU.LEAF)
data$LOG10.CFU.LEAF. <- as.numeric(data$LOG10.CFU.LEAF.)
data$MODA <- paste(data$TREATMENT, data$CONCENTRATION, sep = "_")

Eyach2 <- subset(data,data$ECOTYPE=="Eyach" & data$DPI=="2")
Eyach5 <- subset(data,data$ECOTYPE=="Eyach" & data$DPI=="5")

Col2 <- subset(data,data$ECOTYPE=="Col-0" & data$DPI=="2")
Col5 <- subset(data,data$ECOTYPE=="Col-0" & data$DPI=="5")

desired_order <- c("Mock_-", "WT_10^6", "wfgD_10^6", "WT_10^7", "wfgD_10^7", "WT_10^8", "wfgD_10^8")
Eyach2$MODA <- factor(Eyach2$MODA, levels = desired_order)
Eyach5$MODA <- factor(Eyach5$MODA, levels = desired_order)
Col2$MODA   <- factor(Col2$MODA, levels = desired_order)
Col5$MODA   <- factor(Col5$MODA, levels = desired_order)

# Plot single ecotype - single dpi
fill_colors <- c("grey", "mediumpurple1", "lavender",
                 "mediumpurple1", "lavender", "mediumpurple1", "lavender")
border_colors <- c("black", "purple4", "mediumpurple1",
                   "purple4", "mediumpurple1", "purple4", "mediumpurple1")

make_boxplot <- function(df, title) {
  ggplot(df, aes(x = MODA, y = LOG10.CFU.LEAF.)) +
    geom_boxplot(aes(fill = MODA, color = MODA), width = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = MODA), width = 0.2, size = 1.8, alpha = 0.8) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = border_colors) +
    labs(title = title, x = "Treatment", y = expression(Log[10]~"(CFU/leaf)")) +
    coord_cartesian(ylim = c(0, 8)) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(hjust = 0.5)
    )
}

p1 <- make_boxplot(Eyach2, "Eyach 1.5-2 – 2 dpi")
p2 <- make_boxplot(Eyach5, "Eyach 1.5-2 – 5 dpi")
p3 <- make_boxplot(Col2,   "Col-0 – 2 dpi")
p4 <- make_boxplot(Col5,   "Col-0 – 5 dpi")
p1; p2; p3; p4

# Export plots (PDFs into ./figures)
ggsave("./figures/plant_CFUs_Eyach_2dpi.pdf", p1, width = 6, height = 6)
ggsave("./figures/plant_CFUs_Eyach_5dpi.pdf", p2, width = 6, height = 6)
ggsave("./figures/plant_CFUs_Col0_2dpi.pdf",  p3, width = 6, height = 6)
ggsave("./figures/plant_CFUs_Col0_5dpi.pdf",  p4, width = 6, height = 6)

# Stats
compare_wt_vs_wfgD <- function(df, value_col = "LOG10.CFU.LEAF.") {
  
  results <- data.frame(
    Comparison = character(),
    Test = character(),
    Statistic = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  concentrations <- c("10^6", "10^7", "10^8")
  
  for (conc in concentrations) {
    WT <- paste0("WT_", conc)
    wfgD <- paste0("wfgD_", conc)
    
    # Skip if either group is missing
    if (!all(c(WT, wfgD) %in% df$MODA)) next
    
    x <- df[[value_col]][df$MODA == WT]
    y <- df[[value_col]][df$MODA == wfgD]
    
    # Check normality
    normal_x <- if(length(x) >= 3) shapiro.test(x)$p.value > 0.05 else FALSE
    normal_y <- if(length(y) >= 3) shapiro.test(y)$p.value > 0.05 else FALSE
    
    if (normal_x & normal_y) {
      # Check variance equality
      if (var.test(x, y)$p.value > 0.05) {
        t_res <- t.test(x, y, var.equal = TRUE)
        test_type <- "t-test"
      } else {
        t_res <- t.test(x, y, var.equal = FALSE)
        test_type <- "Welch t-test"
      }
      results <- rbind(results, data.frame(
        Comparison = paste(WT, "vs", wfgD),
        Test = test_type,
        Statistic = t_res$statistic,
        p_value = t_res$p.value
      ))
    } else {
      # Non-normal → Wilcoxon rank-sum
      w_res <- wilcox.test(x, y)
      results <- rbind(results, data.frame(
        Comparison = paste(WT, "vs", wfgD),
        Test = "Wilcoxon rank-sum",
        Statistic = w_res$statistic,
        p_value = w_res$p.value
      ))
    }
  }
  
  return(results)
}

# --- Apply to your datasets ---
res_Eyach2 <- compare_wt_vs_wfgD(Eyach2)
res_Eyach5 <- compare_wt_vs_wfgD(Eyach5)
res_Col2   <- compare_wt_vs_wfgD(Col2)
res_Col5   <- compare_wt_vs_wfgD(Col5)

# Combine results
res_all <- do.call(rbind, list(
  cbind(Dataset = "Eyach2", res_Eyach2),
  cbind(Dataset = "Eyach5", res_Eyach5),
  cbind(Dataset = "Col2",   res_Col2),
  cbind(Dataset = "Col5",   res_Col5)
))

res_all

# Export analysis results
write.csv(res_all, "./data/output/plant_CFUs_stats.csv", row.names = FALSE)
