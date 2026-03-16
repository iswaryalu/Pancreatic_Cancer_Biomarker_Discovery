rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

# =========================
# 1. LOAD LIBRARIES
# =========================
library(glmnet)

# =========================
# 2. LOAD DATA
# =========================
expr  <- read.csv("normalized_expression_matrix.csv", row.names = 1, check.names = FALSE)
genes <- read.csv("Hub_DEG_overlap_genes-feb11.txt")
meta  <- read.csv("meta_clean.csv")

gene_list <- unique(genes$Gene)

# Filter expression matrix
expr_86 <- expr[rownames(expr) %in% gene_list, ]
expr_86_t <- as.data.frame(t(expr_86))

# Sample ID cleanup
expr_86_t$SampleID <- gsub("_.*", "", rownames(expr_86_t))

# Merge metadata
merged_data <- merge(meta, expr_86_t, by = "SampleID")

# Keep Tumor vs Adjacent
merged_data <- merged_data[merged_data$tissue %in% c("Tumor", "adjacent_non_tumor"), ]
merged_data$tissue <- factor(merged_data$tissue)

# Create ML matrix
ml_data <- merged_data[, c("tissue", gene_list)]
ml_data$tissue <- factor(ml_data$tissue)

# Outcome encoding
y <- ifelse(ml_data$tissue == "Tumor", 1, 0)

# Predictor matrix
x <- as.matrix(ml_data[, gene_list])
x <- apply(x, 2, as.numeric)

# Scale properly
x <- scale(x)

# Remove NA samples
complete_idx <- complete.cases(x)
x <- x[complete_idx, ]
y <- y[complete_idx]

cat("Final samples:", nrow(x), "\n")
cat("Final genes:", ncol(x), "\n")

# =========================
# 3. LASSO BIOMARKER DISCOVERY
# =========================
cvfit <- cv.glmnet(
  x, y,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  standardize = FALSE
)

lasso_fit <- glmnet(x, y, family = "binomial", alpha = 1, standardize = FALSE)

cat("Lambda.min:", cvfit$lambda.min, "\n")
cat("Lambda.1se:", cvfit$lambda.1se, "\n")

# Extract genes
coef_lasso <- coef(cvfit, s = "lambda.min")
lasso_genes <- rownames(coef_lasso)[coef_lasso[,1] != 0]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")

cat("Selected biomarker genes:", length(lasso_genes), "\n")
print(lasso_genes)

write.csv(lasso_genes, "~/work/r_project_ishu/ml/lasso/results/LASSO_selected_biomarkers-feb17.csv", row.names = FALSE)

# =========================
# 4. COEFFICIENTS
# =========================
coef_df <- data.frame(
  Gene = lasso_genes,
  Coefficient = coef_lasso[lasso_genes,1]
)
write.csv(coef_df, "~/work/r_project_ishu/ml/lasso/results/LASSO_coefficients-feb17.csv", row.names = FALSE)

# =========================
# 5. PAPER-STYLE LASSO PLOTS
# =========================

dir.create("~/work/r_project_ishu/ml/lasso/plots", recursive = TRUE, showWarnings = FALSE)

# =========================
# 5. IMPROVED PAPER-STYLE LASSO PLOTS
# =========================
library(ggplot2)  # For superior customization [cite:21][cite:22]

dir.create("~/work/r_project_ishu/ml/lasso/plots", recursive = TRUE, showWarnings = FALSE)

# PANEL A: Coefficient Paths (fixed margins, no overcrowding)
png("~/work/r_project_ishu/ml/lasso/plots/LASSO_Coefficient_Path.png", 
    width=2400, height=1800, res=300)  # Larger for labels

fit <- glmnet(x, y, family = "binomial", alpha = 1)

coef_mat <- as.matrix(fit$beta)
log_lambda <- log(fit$lambda)

# Match lengths safely
common_len <- min(ncol(coef_mat), length(log_lambda))
coef_mat <- coef_mat[, 1:common_len]
log_lambda <- log_lambda[1:common_len]
df_vals <- fit$df[1:common_len]   # ⭐ feature count

# Descending
idx <- order(log_lambda, decreasing = TRUE)

matplot(log_lambda[idx],
        t(coef_mat[, idx]),
        type = "l",
        lty = 1,
        lwd = 2,
        xlab = "Log Lambda",
        ylab = "Coefficients",
        main = "LASSO Coefficient Profiles")


# ---- TOP NUMBERS ONLY (like paper) ----
pretty_ticks <- pretty(log_lambda[idx])

# Find closest lambda index for each tick
tick_idx <- sapply(pretty_ticks, function(v) which.min(abs(log_lambda - v)))

# Feature counts
top_numbers <- df_vals[tick_idx]

# Draw numbers only (no label text)
axis(3,
     at = pretty_ticks,
     labels = top_numbers,
     tick = FALSE,
     line = -0.5,   # moves inside box (paper style)
     cex.axis = 1.2)

dev.off()

# PANEL B: CV Deviance (ggplot2 style)
png("~/work/r_project_ishu/ml/lasso/plots/LASSO_Panel_B_CV.png", width=2000, height=1500, res=300)

log_lambda <- log(cvfit$lambda)
idx <- order(log_lambda, decreasing = TRUE)

# Plot CV curve (descending)
plot(log_lambda[idx], cvfit$cvm[idx],
     type = "b",
     pch = 16,
     col = "red",
     xlab = "Log(λ)",
     ylab = "Binomial Deviance",
     main = "LASSO Cross-Validation (Descending Axis)")

# Error bars
arrows(log_lambda[idx],
       cvfit$cvm[idx] - cvfit$cvsd[idx],
       log_lambda[idx],
       cvfit$cvm[idx] + cvfit$cvsd[idx],
       angle = 90, code = 3, length = 0.03, col = "grey60")

# Lambda lines
abline(v = log(cvfit$lambda.min), lty=2)
abline(v = log(cvfit$lambda.1se), lty=2)

# Top numbers (features)
pretty_ticks <- pretty(log_lambda[idx])
tick_idx <- sapply(pretty_ticks, function(v) which.min(abs(log_lambda - v)))
top_numbers <- cvfit$nzero[tick_idx]

axis(3, at = pretty_ticks, labels = top_numbers,
     tick = FALSE, line = -0.5, cex.axis = 1.2)

dev.off()


# PANEL C: Genes vs Log Lambda (ggplot2 for smoothness)
png("~/work/r_project_ishu/ml/lasso/plots/LASSO_Panel_C_Features.png", width=2000, height=1500, res=300)

plot(log_lambda[idx], cvfit$nzero[idx],
     type = "l",
     lwd = 3,
     xlab = "Log(λ)",
     ylab = "Number of Features",
     main = "LASSO Sparsity Profile")

points(log_lambda[idx], cvfit$nzero[idx], pch=16)

abline(v = log(cvfit$lambda.min), lty=2)
abline(v = log(cvfit$lambda.1se), lty=2)

dev.off()

cat("log(lambda.min):", log(cvfit$lambda.min), "\n")
cat("log(lambda.1se):", log(cvfit$lambda.1se), "\n")
