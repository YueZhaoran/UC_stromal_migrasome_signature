##################################################
# Script: 08a_machine_learning_parameters.R
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

# Quick extraction of ML parameters
# Working directory should be configured from project root

library(glmnet)
library(caret)
library(randomForest)
library(e1071)
library(dplyr)

# Load data
dat_bulk <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
group_info <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(group_info) <- c("sample", "group")
hubgene <- read.csv("Venn_intersection_genes.csv", stringsAsFactors = FALSE)
candidate_genes <- hubgene$gene
train_expr <- t(dat_bulk[candidate_genes, ]) %>% as.data.frame()
train_expr$sample <- rownames(train_expr)
train_expr <- merge(train_expr, group_info, by = "sample")
rownames(train_expr) <- train_expr$sample
train_expr <- train_expr[, !(colnames(train_expr) %in% c("sample"))]
train_label <- factor(train_expr$group, levels = c("Normal", "UC"))
train_expr <- train_expr[, !(colnames(train_expr) %in% c("group"))]

# LASSO
set.seed(798)
res.lasso <- cv.glmnet(as.matrix(train_expr), train_label, family = "binomial", type.measure = "default", nfolds = 5)
coef.min <- coef(res.lasso, s = "lambda.min")
lasso_genes <- coef.min@Dimnames[[1]][coef.min@i + 1]
lasso_genes <- lasso_genes[-1]
cat("LASSO lambda.min:", res.lasso$lambda.min, "/n")
cat("LASSO n genes:", length(lasso_genes), "/n")
cat("LASSO genes:", paste(lasso_genes, collapse = ", "), "/n")

# SVM-RFE
set.seed(798)
control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
results <- rfe(train_expr, train_label, sizes = 1:ncol(train_expr), rfeControl = control, method = "svmRadial")
svm_genes <- predictors(results)
cat("SVM-RFE optimal size:", results$optsize, "/n")
cat("SVM-RFE n genes:", length(svm_genes), "/n")
cat("SVM-RFE genes:", paste(svm_genes, collapse = ", "), "/n")

# RF
set.seed(798)
cv_error <- c()
ntree_range <- seq(10, 200, by = 10)
for (ntree in ntree_range) {
  rf_model <- randomForest(x = train_expr, y = train_label, ntree = ntree, mtry = round(sqrt(ncol(train_expr))))
  cv_error <- c(cv_error, rf_model$err.rate[ntree, "OOB"])
}
cv_results <- data.frame(ntree = ntree_range, error = cv_error)
optimal_ntree <- cv_results$ntree[which.min(cv_results$error)]
cat("RF optimal ntree:", optimal_ntree, "/n")

final_rf <- randomForest(x = train_expr, y = train_label, ntree = optimal_ntree, importance = TRUE)
importance_df <- as.data.frame(importance(final_rf))
importance_df$symbol <- rownames(importance_df)
importance_df <- importance_df %>% arrange(desc(MeanDecreaseGini)) %>% filter(MeanDecreaseGini > 0)
cat("RF n genes (>0):", nrow(importance_df), "/n")

# Intersection
lasso_df <- data.frame(symbol = lasso_genes)
svm_df <- data.frame(symbol = svm_genes)
rf_genes <- importance_df$symbol
final_intersection <- Reduce(intersect, list(lasso_df$symbol, rf_genes, svm_df$symbol))
cat("Intersection n genes:", length(final_intersection), "/n")
cat("Intersection genes:", paste(final_intersection, collapse = ", "), "/n")
