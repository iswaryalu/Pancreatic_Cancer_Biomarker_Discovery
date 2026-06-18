rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

# =========================
# 1. LOAD LIBRARIES
# =========================
library(glmnet)
library(ggplot2)

# =========================
# 2. LOAD DATA
# =========================
expr  <- read.csv("~/work/r_project_ishu/deg/data/normalized_expression_matrix.csv",
                  row.names = 1, check.names = FALSE)

genes <- read.delim("~/work/r_project_ishu/wgcna/results/feb11 results/Hub_DEG_overlap_genes-feb11.txt",
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    sep = "\t")

meta  <- read.csv("~/work/r_project_ishu/deg/data/meta_clean.csv")

# =========================
# 3. FIX GENE LIST (IMPORTANT)
# =========================

# Always take first column (your column is 'x')
gene_list <- unique(genes[[1]])

cat("Total genes loaded:", length(gene_list), "\n")
print(head(gene_list))

# =========================
# 4. PREPARE EXPRESSION MATRIX
# =========================

# Filter genes
expr_sub <- expr[rownames(expr) %in% gene_list, ]

# Transpose
expr_t <- as.data.frame(t(expr_sub))

# Add SampleID
expr_t$SampleID <- gsub("_.*", "", rownames(expr_t))

# Merge metadata
merged_data <- merge(meta, expr_t, by = "SampleID")

# Keep Tumor vs Adjacent
merged_data <- merged_data[merged_data$tissue %in% c("Tumor", "adjacent_non_tumor"), ]

# =========================
# 5. SAFE FEATURE SELECTION
# =========================

valid_genes <- intersect(gene_list, colnames(merged_data))

cat("Genes requested:", length(gene_list), "\n")
cat("Genes found:", length(valid_genes), "\n")

# Create ML dataset safely
ml_data <- merged_data[, c("tissue", valid_genes), drop = FALSE]

# Encode outcome
ml_data$tissue <- factor(ml_data$tissue)

y <- ifelse(ml_data$tissue == "Tumor", 1, 0)

# Predictor matrix
x <- as.matrix(ml_data[, valid_genes])

# Convert to numeric safely
x <- apply(x, 2, as.numeric)

# Scale
x <- scale(x)

# Remove NA samples
complete_idx <- complete.cases(x)
x <- x[complete_idx, ]
y <- y[complete_idx]

cat("Final samples:", nrow(x), "\n")
cat("Final genes:", ncol(x), "\n")

# =========================
# 6. LASSO MODEL
# =========================

cvfit <- cv.glmnet(
  x, y,
  family = "binomial",
  alpha = 1,
  nfolds = 10,
  standardize = FALSE
)

fit <- glmnet(x, y, family = "binomial", alpha = 1, standardize = FALSE)

cat("Lambda.min:", cvfit$lambda.min, "\n")
cat("Lambda.1se:", cvfit$lambda.1se, "\n")

# =========================
# 7. EXTRACT GENES
# =========================

coef_lasso <- coef(cvfit, s = "lambda.min")

lasso_genes <- rownames(coef_lasso)[coef_lasso[,1] != 0]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")

cat("Selected genes:", length(lasso_genes), "\n")
print(lasso_genes)

# Save results
write.csv(lasso_genes,
          "~/work/r_project_ishu/ml/lasso/results/LASSO_selected_biomarkers.csv",
          row.names = FALSE)

coef_df <- data.frame(
  Gene = lasso_genes,
  Coefficient = coef_lasso[lasso_genes,1]
)

write.csv(coef_df,
          "~/work/r_project_ishu/ml/lasso/results/LASSO_coefficients.csv",
          row.names = FALSE)

# =========================
# 8. PLOTS
# =========================

dir.create("~/work/r_project_ishu/ml/lasso/plots", recursive = TRUE, showWarnings = FALSE)

# ---- PANEL A: Coefficient Path ----
png("~/work/r_project_ishu/ml/lasso/plots/LASSO_Coefficient_Path.png",
    width = 2400, height = 1800, res = 300)

coef_mat <- as.matrix(fit$beta)
log_lambda <- log(fit$lambda)

common_len <- min(ncol(coef_mat), length(log_lambda))
coef_mat <- coef_mat[, 1:common_len]
log_lambda <- log_lambda[1:common_len]
df_vals <- fit$df[1:common_len]

idx <- order(log_lambda, decreasing = TRUE)

matplot(log_lambda[idx],
        t(coef_mat[, idx]),
        type = "l",
        lty = 1,
        lwd = 2,
        xlab = "Log Lambda",
        ylab = "Coefficients",
        main = "LASSO Coefficient Profiles")

pretty_ticks <- pretty(log_lambda[idx])
tick_idx <- sapply(pretty_ticks, function(v) which.min(abs(log_lambda - v)))
top_numbers <- df_vals[tick_idx]

axis(3, at = pretty_ticks, labels = top_numbers,
     tick = FALSE, line = -0.5, cex.axis = 1.2)

dev.off()

# ---- PANEL B: CV Plot ----
png("~/work/r_project_ishu/ml/lasso/plots/LASSO_CV_Plot.png",
    width = 2000, height = 1500, res = 300)

log_lambda <- log(cvfit$lambda)
idx <- order(log_lambda, decreasing = TRUE)

plot(log_lambda[idx], cvfit$cvm[idx],
     type = "b",
     pch = 16,
     col = "red",
     xlab = "Log(λ)",
     ylab = "Binomial Deviance",
     main = "LASSO Cross-Validation")

arrows(log_lambda[idx],
       cvfit$cvm[idx] - cvfit$cvsd[idx],
       log_lambda[idx],
       cvfit$cvm[idx] + cvfit$cvsd[idx],
       angle = 90, code = 3, length = 0.03, col = "grey60")

abline(v = log(cvfit$lambda.min), lty = 2)
abline(v = log(cvfit$lambda.1se), lty = 2)

pretty_ticks <- pretty(log_lambda[idx])
tick_idx <- sapply(pretty_ticks, function(v) which.min(abs(log_lambda - v)))
top_numbers <- cvfit$nzero[tick_idx]

axis(3, at = pretty_ticks, labels = top_numbers,
     tick = FALSE, line = -0.5, cex.axis = 1.2)

dev.off()

# ---- PANEL C: Feature Count ----
png("~/work/r_project_ishu/ml/lasso/plots/LASSO_Feature_Count.png",
    width = 2000, height = 1500, res = 300)

plot(log_lambda[idx], cvfit$nzero[idx],
     type = "l",
     lwd = 3,
     xlab = "Log(λ)",
     ylab = "Number of Features",
     main = "LASSO Sparsity Profile")

points(log_lambda[idx], cvfit$nzero[idx], pch = 16)

abline(v = log(cvfit$lambda.min), lty = 2)
abline(v = log(cvfit$lambda.1se), lty = 2)

dev.off()

# =========================
# 9. FINAL OUTPUT
# =========================
cat("log(lambda.min):", log(cvfit$lambda.min), "\n")
cat("log(lambda.1se):", log(cvfit$lambda.1se), "\n")