##################################################
# Script: 13b_network_statistics.R
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

# Quick extraction of ANN performance
# Working directory should be configured from project root

library(neuralnet)
library(pROC)
library(ggplot2)

hub_genes <- read.csv("intersection_gene.csv", stringsAsFactors = FALSE)
key_genes <- hub_genes$symbol

train_data <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
train_group <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(train_group) <- c("sample", "group")

val_data <- read.csv("../00bulkdata/dat.GSE75214.csv", row.names = 1, check.names = FALSE)
val_group <- read.csv("../00bulkdata/group.GSE75214.csv", stringsAsFactors = FALSE)
colnames(val_group) <- c("sample", "group")

train_df <- as.data.frame(t(train_data[key_genes, ]))
train_df$group <- train_group$group
train_df$label <- ifelse(train_df$group == "UC", 1, 0)

val_df <- as.data.frame(t(val_data[key_genes, ]))
val_df$group <- val_group$group
val_df$label <- ifelse(val_df$group == "UC", 1, 0)

train_df_scaled <- train_df
train_df_scaled[, key_genes] <- scale(train_df_scaled[, key_genes])

formula_str <- paste("label ~", paste(key_genes, collapse = " + "))
formula_obj <- as.formula(formula_str)

set.seed(798)
ann_model <- neuralnet(formula_obj, data = train_df_scaled, hidden = 5, linear.output = FALSE)

train_pred <- compute(ann_model, train_df_scaled[, key_genes])
roc_train <- roc(train_df_scaled$label, train_pred$net.result[, 1], quiet = TRUE)
cat("Training AUC:", round(auc(roc_train), 3), "/n")

val_df_scaled <- val_df
val_df_scaled[, key_genes] <- scale(val_df_scaled[, key_genes])
val_pred <- compute(ann_model, val_df_scaled[, key_genes])
roc_val <- roc(val_df_scaled$label, val_pred$net.result[, 1], quiet = TRUE)
cat("Validation AUC:", round(auc(roc_val), 3), "/n")

# Additional metrics at optimal threshold
get_metrics <- function(true_labels, pred_probs) {
  coords_res <- coords(roc(true_labels, pred_probs, quiet = TRUE), "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  threshold <- coords_res$threshold[1]
  pred_labels <- ifelse(pred_probs >= threshold, 1, 0)
  cm <- table(Predicted = pred_labels, Actual = true_labels)
  accuracy <- sum(diag(cm)) / sum(cm)
  cat("Threshold:", round(threshold, 3), "/n")
  cat("Sensitivity:", round(coords_res$sensitivity[1], 3), "/n")
  cat("Specificity:", round(coords_res$specificity[1], 3), "/n")
  cat("Accuracy:", round(accuracy, 3), "/n")
}

cat("/n=== Training metrics ===/n")
get_metrics(train_df_scaled$label, train_pred$net.result[, 1])

cat("/n=== Validation metrics ===/n")
get_metrics(val_df_scaled$label, val_pred$net.result[, 1])
