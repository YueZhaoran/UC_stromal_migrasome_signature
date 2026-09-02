##################################################
# Script: 09_model_validation.R
#
# Purpose:
# Analysis workflow for the study:
# Integrating genetic risk and single-cell transcriptomics
# identifies a stromal migrasome-related signature in ulcerative colitis
#
# This script was prepared for reproducible research.
# Please configure input/output paths before execution.
##################################################
# Input files should be placed under data/ and resources/.
# Output files should be written under results/.
##################################################

# Quick extraction of verification stats
# Working directory should be configured from project root

library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggh4x)
library(pROC)
library(tidyr)
library(dplyr)

hub_gene <- read.csv("intersection_gene.csv", stringsAsFactors = FALSE)

# Training
train_data <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
train_group <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(train_group) <- c("sample", "group")

temp1 <- train_data[rownames(train_data) %in% hub_gene$symbol, ]
temp1$symbol <- rownames(temp1)
temp1_long <- pivot_longer(temp1, cols = -symbol, names_to = 'sample', values_to = 'expr')
train_hub_exp <- merge(temp1_long, train_group, by = 'sample')
train_hub_exp$group <- as.factor(train_hub_exp$group)

train_stat.test <- train_hub_exp %>%
  group_by(symbol) %>%
  wilcox_test(expr ~ group) %>%
  adjust_pvalue(method = 'fdr')

cat("=== Training set (GSE87466) ===/n")
print(train_stat.test[, c("symbol", "p", "p.adj")])

biomarker_exp_train <- t(train_data[hub_gene$symbol, ])
for (gene in hub_gene$symbol) {
  roc_obj <- roc(train_group$group, biomarker_exp_train[, gene], levels = c('Normal', 'UC'), quiet = TRUE)
  cat(gene, "AUC:", round(roc_obj$auc, 3), "/n")
}

# Validation
ver_data <- read.csv("../00bulkdata/dat.GSE75214.csv", row.names = 1, check.names = FALSE)
ver_group <- read.csv("../00bulkdata/group.GSE75214.csv", stringsAsFactors = FALSE)
colnames(ver_group) <- c("sample", "group")

temp1_ver <- ver_data[rownames(ver_data) %in% hub_gene$symbol, ]
temp1_ver$symbol <- rownames(temp1_ver)
temp1_ver_long <- pivot_longer(temp1_ver, cols = -symbol, names_to = 'sample', values_to = 'expr')
ver_hub_exp <- merge(temp1_ver_long, ver_group, by = 'sample')
ver_hub_exp$group <- as.factor(ver_hub_exp$group)

ver_stat.test <- ver_hub_exp %>%
  group_by(symbol) %>%
  wilcox_test(expr ~ group) %>%
  adjust_pvalue(method = 'fdr')

cat("/n=== Validation set (GSE75214) ===/n")
print(ver_stat.test[, c("symbol", "p", "p.adj")])

biomarker_exp_ver <- t(ver_data[hub_gene$symbol, ])
for (gene in hub_gene$symbol) {
  roc_obj <- roc(ver_group$group, biomarker_exp_ver[, gene], levels = c('Normal', 'UC'), quiet = TRUE)
  cat(gene, "AUC:", round(roc_obj$auc, 3), "/n")
}

# Consistency check
disease_sample <- train_group[train_group$group == 'UC', ]
control_sample <- train_group[train_group$group == 'Normal', ]
train_mean <- data.frame(
  disease_mean = rowMeans(temp1[, disease_sample$sample]),
  control_mean = rowMeans(temp1[, control_sample$sample])
)
train_mean$compare_res <- ifelse(train_mean$disease_mean > train_mean$control_mean, 1, 0)
train_mean$symbol <- rownames(train_mean)

disease_sample_ver <- ver_group[ver_group$group == 'UC', ]
control_sample_ver <- ver_group[ver_group$group == 'Normal', ]
ver_mean <- data.frame(
  disease_mean = rowMeans(temp1_ver[, disease_sample_ver$sample]),
  control_mean = rowMeans(temp1_ver[, control_sample_ver$sample])
)
ver_mean$compare_res <- ifelse(ver_mean$disease_mean > ver_mean$control_mean, 1, 0)
ver_mean$symbol <- rownames(ver_mean)

cat("/n=== Expression trend (1=UC high, 0=UC low) ===/n")
print(merge(train_mean[, c("symbol", "compare_res")], ver_mean[, c("symbol", "compare_res")], by = "symbol"))
