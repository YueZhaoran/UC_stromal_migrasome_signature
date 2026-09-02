##################################################
# Script: 08_machine_learning_signature.R
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
# 设置项目目录
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
project_dir  <- dirname(rstudioapi::getActiveDocumentContext()$path)
cat("工作目录：", getwd(), "/n")


library(glmnet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomForest)
library(caret)
library(e1071)
library(ggvenn)


# 加载数据与候选基因


# 读取 Bulk 表达矩阵
dat_bulk <- read.csv("../00bulkdata/dat.GSE87466.csv", row.names = 1, check.names = FALSE)

# 读取分组信息
group_info <- read.csv("../00bulkdata/group.GSE87466.csv", stringsAsFactors = FALSE)
colnames(group_info) <- c("sample", "group")  # 统一列名

# 读取韦恩图取到的候选基因
hubgene <- read.csv("Venn_intersection_genes.csv", stringsAsFactors = FALSE)
candidate_genes <- hubgene$gene

# 提取表达矩阵并转置
train_expr <- t(dat_bulk[candidate_genes, ]) %>% as.data.frame()

# 确保分组标签与表达矩阵的样本顺序一致
train_expr$sample <- rownames(train_expr)
train_expr <- merge(train_expr, group_info, by = "sample")
rownames(train_expr) <- train_expr$sample
train_expr <- train_expr %>% dplyr::select(-sample)

# 将分组转为 0/1 因子（1=UC, 0=Normal）
train_label <- factor(train_expr$group, levels = c("Normal", "UC"))
train_expr <- train_expr %>% dplyr::select(-group)

cat(paste0("✅ 成功加载 ", nrow(train_expr), " 个样本, ", ncol(train_expr), " 个候选基因。/n"))

# 创建结果输出文件夹



# ==================================================================
# LASSO 回归筛选
# ==================================================================

set.seed(798)
res.lasso <- cv.glmnet(as.matrix(train_expr), train_label, 
                       family = "binomial", type.measure = "default", nfolds = 5)
coef.min <- coef(res.lasso, s = "lambda.min")
active.min <- which(coef.min@i != 0)
lasso_genes <- coef.min@Dimnames[[1]][coef.min@i + 1]
lasso_genes <- lasso_genes[-1]  # 去掉截距项
write.csv(data.frame(symbol = lasso_genes), "01.lasso_gene.csv", row.names = FALSE)

# 绘图1：LASSO 系数路径
fit <- glmnet(as.matrix(train_expr), train_label, family = "binomial")
tmp <- as.data.frame(as.matrix(coef(fit)))
tmp$coef <- rownames(tmp)
tmp <- reshape2::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$lambda <- fit$lambda[tmp$variable + 1]
tmp <- tmp[tmp$coef != "(Intercept)", ]

p_lasso <- ggplot(tmp, aes(log(lambda), value, color = coef)) +
  geom_vline(xintercept = log(res.lasso$lambda.min), linetype = 2, color = "grey60") +
  geom_line(size = 0.8) +
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) +
  labs(x = "Log(lambda)", y = "Coefficients") +
  annotate("text", x = -4, y = -3, label = paste0("Optimal Lambda = ", round(res.lasso$lambda.min, 4)))

ggsave("02.Lasso_model.pdf", p_lasso, width = 8, height = 6)
ggsave("02.Lasso_model.png", p_lasso, width = 8, height = 6)

# 绘图2：LASSO 交叉验证
xx <- data.frame(lambda = res.lasso[["lambda"]], cvm = res.lasso[["cvm"]], 
                 cvup = res.lasso[["cvup"]], cvlo = res.lasso[["cvlo"]], 
                 nozezo = res.lasso[["nzero"]])
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo, " vars")

p_cv <- ggplot(xx, aes(ll, cvm, color = NZERO)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.05) +
  geom_vline(xintercept = xx$ll[which.min(xx$cvm)], linetype = 2, color = "grey60") +
  geom_point(size = 2) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Log(lambda)", y = "Partial Likelihood Deviance") +
  annotate("text", x = -4.5, y = max(xx$cvm) - 0.2, label = paste0("Optimal Lambda = ", round(res.lasso$lambda.min, 4)))
ggsave("03.Lasso_verify.pdf", p_cv, width = 8, height = 6)
ggsave("03.Lasso_verify.png", p_cv, width = 8, height = 6)


# ==================================================================
# SVM-RFE 递归特征消除
# ==================================================================

set.seed(798)
control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
num <- ncol(train_expr)

results <- rfe(train_expr, train_label, sizes = 1:num, rfeControl = control, method = "svmRadial")

svm_genes <- predictors(results)
write.csv(data.frame(symbol = svm_genes), "svm_rfe_gene.csv", row.names = FALSE)

# SVM-RFE 绘图
svm_plot_data <- results$results
svm_plot_data$col <- ifelse(svm_plot_data$Kappa == max(svm_plot_data$Kappa), "1", "0")

p_acc <- ggplot(svm_plot_data, aes(Variables, Accuracy)) +
  geom_line(color = "#4865A9") +
  geom_point(aes(color = col), size = 3) +
  scale_color_manual(values = c("#4865A9", "#EF8A43")) +
  annotate("text", x = svm_plot_data$Variables[svm_plot_data$col == "1"] + 0.5,
           y = max(svm_plot_data$Accuracy) - 0.001, color = "red",
           label = paste0("(", svm_plot_data$Variables[svm_plot_data$col == "1"], ")")) +
  theme_bw() + theme(legend.position = "none")

p_err <- ggplot(svm_plot_data, aes(Variables, 1 - Accuracy)) +
  geom_line(color = "#4865A9") +
  geom_point(aes(color = col), size = 3) +
  scale_color_manual(values = c("#4865A9", "#EF8A43")) +
  theme_bw() + theme(legend.position = "none")

library(patchwork)
p_svm <- p_acc + p_err
ggsave("01.SVM_RFE.pdf", p_svm, width = 14, height = 5)
ggsave("01.SVM_RFE.png", p_svm, width = 14, height = 5)


# ==================================================================
# 随机森林 (Random Forest)
# ==================================================================

set.seed(798)
cv_error <- c()
ntree_range <- seq(10, 200, by = 10)

for (ntree in ntree_range) {
  rf_model <- randomForest(x = train_expr, y = train_label, ntree = ntree, 
                           mtry = round(sqrt(ncol(train_expr))))
  cv_error <- c(cv_error, rf_model$err.rate[ntree, "OOB"])
}
cv_results <- data.frame(ntree = ntree_range, error = cv_error)
optimal_ntree <- cv_results$ntree[which.min(cv_results$error)]

# RF 调优图
p_rf_err <- ggplot(cv_results, aes(ntree, error)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = optimal_ntree, linetype = 2, color = "red") +
  labs(x = "Number of Trees", y = "OOB Error Rate") + theme_bw()
ggsave("01_RF_select_ntree.pdf", p_rf_err, width = 6, height = 4)
ggsave("01_RF_select_ntree.png", p_rf_err, width = 6, height = 4)

# RF 最优模型
final_rf <- randomForest(x = train_expr, y = train_label, ntree = optimal_ntree, importance = TRUE)
importance_df <- as.data.frame(importance(final_rf))
importance_df$symbol <- rownames(importance_df)
importance_df <- importance_df %>% arrange(desc(MeanDecreaseGini)) %>% filter(MeanDecreaseGini > 0)
write.csv(importance_df,"02.RF_importance.csv", row.names = FALSE)

# RF 重要性图（Top 10）
importance_top10 <- head(importance_df, 10)
p_top <- ggplot(importance_top10, aes(x = reorder(symbol, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#DC0000B2", width = 0.7) + coord_flip() +
  labs(x = "Gene", y = "Importance", title = "Feature Importance") + theme_bw()
ggsave("03_RF_importance.pdf", p_top, width = 8, height = 6)
ggsave("03_RF_importance.png", p_top, width = 8, height = 6)


# ==================================================================
# 三种算法取交集并画图
# ==================================================================


# 读取三份结果
lasso_df <- read.csv("01.lasso_gene.csv")
svm_df <- read.csv("svm_rfe_gene.csv")
rf_genes <- importance_df$symbol  # 使用全部重要性>0的基因，而非仅前10

# 求交集
final_intersection <- Reduce(intersect, list(lasso_df$symbol, rf_genes, svm_df$symbol))
write.csv(data.frame(symbol = final_intersection), "intersection_gene.csv", row.names = FALSE)

cat(paste0("✅ 三种机器学习最终锁定 ", length(final_intersection), " 个共享关键基因。/n"))

# 绘制韦恩图
ml_datalist <- list("LASSO" = lasso_df$symbol, "RF" = rf_genes, "SVM_RFE" = svm_df$symbol)
p_ml <- ggvenn(ml_datalist, 
               columns = c("LASSO", "RF", "SVM_RFE"),
               fill_color = c("#ffb2b2", "turquoise", "#D6E7A3"),
               show_percentage = TRUE, 
               text_size = 4)

ggsave("ML_hub_gene.pdf", p_ml, width = 8, height = 8)
ggsave("ML_hub_gene.png", p_ml, width = 8, height = 8, dpi = 600)

