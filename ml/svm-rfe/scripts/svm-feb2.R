library(caret)
library(e1071)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(doParallel)

stopImplicitCluster()
registerDoSEQ()

set.seed(123)

# =============================
# LOAD DATA
# =============================
expr  <- read.csv("normalized_expression_matrix.csv", row.names = 1, check.names = FALSE)
genes <- read.csv("Hub_DEG_overlap_genes-feb11.txt")
meta  <- read.csv("meta_clean.csv", stringsAsFactors = FALSE)

# =============================
# FILTER GENES
# =============================
gene_list <- unique(genes$Gene)
expr_filt <- expr[rownames(expr) %in% gene_list, ]

expr_t <- as.data.frame(t(expr_filt))
expr_t$SampleID <- gsub("_.*", "", rownames(expr_t))

merged <- merge(meta, expr_t, by="SampleID")

merged <- merged[merged$tissue %in% c("Tumor","adjacent_non_tumor"), ]
merged$tissue <- factor(merged$tissue)

ml_data <- merged[, c("tissue", gene_list)]

# Encode labels
y <- factor(ml_data$tissue)

# IMPORTANT: x must be DATA FRAME
x <- as.data.frame(apply(ml_data[, gene_list], 2, as.numeric))
x <- as.data.frame(scale(x))

# Remove NA
idx <- complete.cases(x)
x <- x[idx,]
y <- y[idx]

cat("Samples:", nrow(x), " Genes:", ncol(x), "\n")

# =========================
# SVM-RFE
# =========================
ctrl <- rfeControl(functions = caretFuncs,
                   method = "cv",
                   number = 5)

sizes <- seq(1, min(50, ncol(x)), by = 1)

svmProfile <- rfe(x, y,
                  sizes = sizes,
                  rfeControl = ctrl,
                  method = "svmRadial",
                  tuneLength = 5)

# =========================
# Extract results
# =========================
results <- svmProfile$results
results$CV_Error <- 1 - results$Accuracy

optimal_gene_number <- results$Variables[which.max(results$Accuracy)]
cat("Optimal gene number =", optimal_gene_number, "\n")


# Optimal indices
opt_acc_idx <- which.max(results$Accuracy)
opt_err_idx <- which.min(results$CV_Error)

opt_features_acc <- results$Variables[opt_acc_idx]
opt_features_err <- results$Variables[opt_err_idx]

opt_acc <- round(results$Accuracy[opt_acc_idx], 3)
opt_err <- round(results$CV_Error[opt_err_idx], 4)

cat("Optimal features (Accuracy) =", opt_features_acc, "\n")
cat("Optimal features (Error) =", opt_features_err, "\n")

# =========================
# 📈 1️⃣ ACCURACY PLOT
# =========================
png("~/work/r_project_ishu/ml/svm-rfe/plots/SVM_RFE_Accuracy.png",
    width = 1600, height = 1200, res = 200)

par(mar = c(5,5,4,2))

# Dynamic headroom
acc_range <- range(results$Accuracy)
acc_pad <- diff(acc_range) * 0.08

plot(results$Variables,
     results$Accuracy,
     type = "l",
     col = "#4C84C4",
     lwd = 3,
     ylim = c(acc_range[1], acc_range[2] + acc_pad),
     xlab = "Number of features",
     ylab = "Accuracy")

points(opt_features_acc, opt_acc, pch = 19)

text(opt_features_acc,
     opt_acc + acc_pad * 0.5,
     labels = paste0(opt_features_acc, " – ", opt_acc),
     col = "red",
     cex = 1.3)

dev.off()

# =========================
# 📉 2️⃣ ERROR PLOT
# =========================
png("~/work/r_project_ishu/ml/svm-rfe/plots/SVM_RFE_Error.png",
    width = 1600, height = 1200, res = 200)

par(mar = c(5,5,4,2))

err_range <- range(results$CV_Error)
err_pad <- diff(err_range) * 0.08

plot(results$Variables,
     results$CV_Error,
     type = "l",
     col = "#4C84C4",
     lwd = 3,
     ylim = c(err_range[1] - err_pad, err_range[2] + err_pad),
     xlab = "Number of features",
     ylab = "5× CV error")

points(opt_features_err, opt_err, pch = 19)

text(opt_features_err,
     opt_err - err_pad * 0.5,
     labels = paste0(opt_features_err, " – ", opt_err),
     col = "red",
     cex = 1.3)

dev.off()


# =========================
# Extract Biomarker Genes
# =========================
svm_genes <- svmProfile$optVariables
write.csv(svm_genes,
          "~/work/r_project_ishu/ml/svm-rfe/results/SVM_RFE_biomarker_genes-feb18.csv",
          row.names = FALSE)
