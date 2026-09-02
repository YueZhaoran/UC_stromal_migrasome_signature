##################################################
# Script: 16_additional_analysis_VNN.R
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

rm(list = ls())
gc()

# ============================ 环境准备 ============================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir  <- dirname(rstudioapi::getActiveDocumentContext()$path)
cat("工作目录：", getwd(), "/n")



# 加载需要使用的包
library(neuralnet)
library(pROC)
library(ggplot2)

# 读取数据 
# 读取最终机器学习选出的关键基因
hub_genes <- read.csv("intersection_gene.csv", stringsAsFactors = FALSE)
key_genes <- hub_genes$symbol

# 读取训练集数据
train_data <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)
train_group <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(train_group) <- c("sample", "group")

# 读取验证集数据
val_data <- read.csv("../00bulkdata/dat.GSE75214.csv", row.names = 1, check.names = FALSE)
val_group <- read.csv("../00bulkdata/group.GSE75214.csv", stringsAsFactors = FALSE)
colnames(val_group) <- c("sample", "group")

# 提取这 3 个基因的表达矩阵，并转置为 行=样本，列=基因
train_df <- as.data.frame(t(train_data[key_genes, ]))
train_df$group <- train_group$group
train_df$label <- ifelse(train_df$group == "UC", 1, 0) 

val_df <- as.data.frame(t(val_data[key_genes, ]))
val_df$group <- val_group$group
val_df$label <- ifelse(val_df$group == "UC", 1, 0)


# 构建人工神经网络（ANN）

train_df_scaled <- train_df
train_df_scaled[, key_genes] <- scale(train_df_scaled[, key_genes])

# 生成模型公式：Label ~ 基因1 + 基因2 + 基因3
formula_str <- paste("label ~", paste(key_genes, collapse = " + "))
formula_obj <- as.formula(formula_str)

# 训练 ANN 模型（设置一个隐藏层，节点数为 5）
set.seed(798)
ann_model <- neuralnet(formula_obj, 
                       data = train_df_scaled, 
                       hidden = 5, 
                       linear.output = FALSE)

# 绘制并保存网络拓扑结构图（图5A）
pdf("ANN_Network_Structure.pdf", width = 10, height = 8)
plot(ann_model, rep = "best", intercept = FALSE, information = FALSE)
dev.off()
png("ANN_Network_Structure.png", width = 10, height = 8, units = "in", res = 600)
plot(ann_model, rep = "best", intercept = FALSE, information = FALSE)
dev.off()


# ============================ 3. 基于 ANN 模型做 ROC 预测 ============================
# --- 训练集预测与 ROC ---
train_pred <- compute(ann_model, train_df_scaled[, key_genes])
roc_train <- roc(train_df_scaled$label, train_pred$net.result[, 1], quiet = TRUE)

pdf("Training_ANN_ROC.pdf", width = 6, height = 6)
plot(roc_train, col = "#D55E00", main = "Training Set ANN ROC", lwd = 2)
legend("bottomright", legend = paste0("AUC = ", round(auc(roc_train), 3)), col = "#D55E00", lwd = 2)
dev.off()
png("Training_ANN_ROC.png", width = 6, height = 6, units = "in", res = 600)
plot(roc_train, col = "#D55E00", main = "Training Set ANN ROC", lwd = 2)
legend("bottomright", legend = paste0("AUC = ", round(auc(roc_train), 3)), col = "#D55E00", lwd = 2)
dev.off()


# --- 验证集预测与 ROC ---
val_df_scaled <- val_df
val_df_scaled[, key_genes] <- scale(val_df_scaled[, key_genes])

val_pred <- compute(ann_model, val_df_scaled[, key_genes])
roc_val <- roc(val_df_scaled$label, val_pred$net.result[, 1], quiet = TRUE)

pdf("Validation_ANN_ROC.pdf", width = 6, height = 6)
plot(roc_val, col = "#0072B2", main = "Validation Set ANN ROC", lwd = 2)
legend("bottomright", legend = paste0("AUC = ", round(auc(roc_val), 3)), col = "#0072B2", lwd = 2)
dev.off()
png("Validation_ANN_ROC.png", width = 6, height = 6, units = "in", res = 600)
plot(roc_val, col = "#0072B2", main = "Validation Set ANN ROC", lwd = 2)
legend("bottomright", legend = paste0("AUC = ", round(auc(roc_val), 3)), col = "#0072B2", lwd = 2)
dev.off()
