##################################################
# Script: 12b_immune_statistics.R
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

# Quick extraction of immune microenvironment stats
# Working directory should be configured from project root

library(dplyr)
library(tidyr)

xcell_score <- read.csv("01.xCell_scores.csv", row.names = 1, check.names = FALSE)
train_data <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
hub_gene <- read.csv("../05DEGs/intersection_gene.csv", stringsAsFactors = FALSE)
key_genes <- hub_gene$symbol

# Direction of DE cells
cell_means <- xcell_score %>%
  group_by(group) %>%
  summarise(across(-any_of("group"), ~mean(.x, na.rm = TRUE))) %>%
  t() %>% as.data.frame()
cell_means <- cell_means[-1, ]
colnames(cell_means) <- c("Normal", "UC")
cell_means$Normal <- as.numeric(cell_means$Normal)
cell_means$UC <- as.numeric(cell_means$UC)
cell_means$Cell <- rownames(cell_means)
cell_means$direction <- ifelse(cell_means$UC > cell_means$Normal, "UC high", "Normal high")

wilcox_res <- read.csv("02.wilcox_res.csv", stringsAsFactors = FALSE)
DE_cells <- read.csv("03.DE_cells.csv", stringsAsFactors = FALSE)

cat("Total DE cells:", nrow(DE_cells), "/n")
cat("/nTop 10 significant DE cells:/n")
DE_cells_sorted <- DE_cells[order(DE_cells$p.adj), ]
for (i in 1:10) {
  cell <- DE_cells_sorted$ImmuneCell[i]
  dir <- cell_means$direction[cell_means$Cell == cell]
  cat(sprintf("%s: padj=%.2e, direction=%s/n", cell, DE_cells_sorted$p.adj[i], dir))
}

# Correlation stats
cor_res <- read.csv("04.Gene_Cell_Correlation.csv", stringsAsFactors = FALSE)
cat("/nTotal gene-cell pairs:", nrow(cor_res), "/n")
cat("Significant pairs (P < 0.05):", sum(cor_res$P < 0.05), "/n")

# Top correlations
sig_cor <- cor_res[cor_res$P < 0.05, ]
sig_cor <- sig_cor[order(abs(sig_cor$Cor), decreasing = TRUE), ]
cat("/nTop 10 strongest correlations:/n")
for (i in 1:10) {
  cat(sprintf("%s - %s: cor=%.3f, p=%.2e/n", sig_cor$Gene[i], sig_cor$Cell[i], sig_cor$Cor[i], sig_cor$P[i]))
}

# Per-gene significant counts
cat("/nSignificant correlations per gene:/n")
for (gene in key_genes) {
  n_sig <- sum(cor_res$P < 0.05 & cor_res$Gene == gene)
  cat(gene, ":", n_sig, "/n")
}
