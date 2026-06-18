# ======================================================
# GSE15932 – External Validation (Blood)
# 4-Gene Logistic Regression Panel
# ======================================================

rm(list=ls())
graphics.off()
set.seed(123)

# =========================
# 1️⃣ Load Libraries
# =========================
library(pROC)
library(ggpubr)

# =========================
# 2️⃣ Define Biomarkers
# =========================
final_genes <- c("AHNAK2","TRIM29","SLC6A14","ITGB4","CTSE")

# =========================
# 3️⃣ Load Data
# =========================
expr <- read.csv(
  "~/work/r_project_ishu/validation_blood/dataset-2/data/normalized_expression_matrix_GSE15932.csv",
  row.names = 1,
  check.names = FALSE
)

meta <- read.csv(
  "~/work/r_project_ishu/validation_blood/dataset-2/data/metadata_GSE15932_clean.csv"
)

# =========================
# 4️⃣ Recode Groups
# =========================
# =========================
# 4️⃣ Subset ONLY PC vs Healthy
# =========================
meta$group <- NA

meta$group[meta$Group == "PC only"] <- "Tumor"
meta$group[meta$Group == "Healthy"] <- "Control"

# Remove all other samples (PC+Diabetes, Diabetes only)
meta <- meta[!is.na(meta$group), ]

meta$group <- factor(meta$group, levels=c("Control","Tumor"))

cat("Filtered Group distribution:\n")
print(table(meta$group))

# =========================
# 5️⃣ Extract Biomarker Expression Safely
# =========================
common_genes <- intersect(final_genes, rownames(expr))
cat("Genes used:", common_genes, "\n")

expr_t <- as.data.frame(t(expr[common_genes, ]))
expr_t$SampleID <- rownames(expr_t)

# Clean Sample IDs
meta$SampleID  <- trimws(as.character(meta$SampleID))
expr_t$SampleID <- trimws(as.character(expr_t$SampleID))

# Merge
df <- merge(meta, expr_t, by="SampleID")

df$group <- factor(df$group, levels=c("Control","Tumor"))

cat("Final dataset size:", nrow(df), "samples\n")
print(table(df$group))

# =========================
# 6️⃣ Logistic Regression Model
# =========================
glm_model <- glm(group ~ .,
                 data = df[, c("group", common_genes)],
                 family = binomial(link="logit"))

cat("\nModel Summary:\n")
print(summary(glm_model))

# =========================
# 7️⃣ Predicted Risk Scores
# =========================
pred <- predict(glm_model, type="response")

cat("\nPrediction Summary:\n")
print(summary(pred))

# =========================
# 8️⃣ Multi-Gene ROC Analysis
# =========================
roc_panel <- pROC::roc(response = df$group,
                       predictor = pred)

cat("\nMulti-gene AUC:\n")
print(auc(roc_panel))

cat("\nDeLong CI:\n")
print(ci.auc(roc_panel))

cat("\nBootstrap CI (2000 resamples):\n")
print(ci.auc(roc_panel, method="bootstrap", boot.n=2000))
# =========================
# 8B️⃣ Single-Gene ROC Analysis
# =========================

single_auc <- data.frame(Gene=character(), AUC=numeric())

dir.create("~/work/r_project_ishu/validation_blood/dataset-2/plots/single_gene_ROC_april2",
           showWarnings = FALSE)

for (gene in common_genes) {
  
  roc_single <- roc(df$group, df[[gene]])
  
  auc_val <- auc(roc_single)
  
  single_auc <- rbind(single_auc,
                      data.frame(Gene=gene, AUC=as.numeric(auc_val)))
  
  # Save individual ROC
  png(paste0("~/work/r_project_ishu/validation_blood/dataset-2/plots/single_gene_ROC_april2/",
             gene, "_ROC_april2.png"),
      width=700, height=600)
  
  plot(roc_single,
       col="red",
       lwd=2,
       main=paste(gene, "ROC"))
  
  text(0.6, 0.2,
       labels=paste("AUC =", round(auc_val,3)),
       cex=1.1)
  
  dev.off()
}

# Sort by best gene
single_auc <- single_auc[order(-single_auc$AUC), ]

print(single_auc)
write.csv(single_auc,
          "~/work/r_project_ishu/validation_blood/dataset-2/results/GSE15932_single_gene_AUC_april2.csv",
          row.names=FALSE)

# =========================
# 9️⃣ Plot ROC
# =========================
# =========================
# 9️⃣ Save Multi-Gene ROC (Base)
# =========================
# =========================
# 9️⃣ Save Multi-Gene ROC (ENHANCED)
# =========================
png("~/work/r_project_ishu/validation_blood/dataset-2/plots/GSE15932_multi_gene_ROC_april2.png",
    width=900, height=750)   # 🔥 bigger canvas

# Plot ROC with stronger emphasis
plot(roc_panel,
     col="blue",
     lwd=5,                       # 🔥 thicker curve
     main="GSE15932 Multi-Gene Logistic ROC",
     cex.main=1.5,                # 🔥 bigger title
     
     print.auc=TRUE,
     print.auc.cex=1.5,           # 🔥 bigger AUC text
     print.auc.y=0.25)

# Add diagonal reference line
abline(a=0, b=1, lty=2, col="gray")

# Add confidence interval
ci_val <- ci.auc(roc_panel)

text(0.6, 0.1,
     labels=paste("95% CI:",
                  round(ci_val[1],3), "-",
                  round(ci_val[3],3)),
     cex=1.2)                    # 🔥 bigger CI text

dev.off()
# Overlay ROC curves
# =========================
# =========================
# Overlay ROC curves (IMPROVED + BIGGER LEGEND)
# =========================
png("~/work/r_project_ishu/validation_blood/dataset-2/plots/GSE15932_ROC_overlay_april2.png",
    width=900, height=750)

# ---- Multi-gene ROC (make dominant) ----
plot(roc_panel,
     col="black",
     lwd=5,                          # 🔥 thicker line
     main="Multi vs Single Gene ROC",
     print.auc=TRUE,
     print.auc.cex=1.4,
     print.auc.y=0.25)

# Prepare legend labels with AUC
legend_labels <- c(
  paste0("Multi-gene (AUC = ", round(auc(roc_panel), 3), ")")
)

# Define colors
cols <- c("red", "blue", "green", "purple", "orange")

i <- 1
for (gene in common_genes) {
  
  roc_single <- roc(df$group, df[[gene]])
  auc_val <- auc(roc_single)
  
  # Add single gene ROC (thinner)
  plot(roc_single, add=TRUE, col=cols[i], lwd=2)
  
  # Add AUC to legend
  legend_labels <- c(
    legend_labels,
    paste0(gene, " (AUC = ", round(auc_val, 3), ")")
  )
  
  i <- i + 1
}

# ---- BIGGER LEGEND ----
legend("bottomright",
       legend = legend_labels,
       col = c("black", cols),
       lwd = c(5, rep(2, length(common_genes))),
       
       cex = 1.4,        # 🔥 bigger text
       pt.cex = 1.6,     # 🔥 bigger line symbols
       y.intersp = 1.5,  # 🔥 more vertical spacing
       x.intersp = 1.3,  # 🔥 more horizontal spacing
       seg.len = 3.5,    # 🔥 longer legend lines
       
       box.lwd = 2,      
       inset = 0.02
)

dev.off()

# =========================
# 🔟 Optimal Cutoff (Youden Index)
# =========================
opt_cut <- coords(roc_panel,
                  "best",
                  ret=c("threshold","sensitivity","specificity"),
                  best.method="youden")

cat("\nOptimal Cutoff (Youden):\n")
print(opt_cut)

# =========================
# 11️⃣ Confusion Matrix at Optimal Cutoff
# =========================
# Select one threshold
threshold <- opt_cut$threshold[1]   # choose row 1

pred_class <- ifelse(pred >= threshold, "Tumor", "Control")
pred_class <- factor(pred_class, levels=c("Control","Tumor"))

cat("\nConfusion Matrix:\n")
print(table(Predicted=pred_class, Actual=df$group))

cat("\nConfusion Matrix:\n")
print(table(Predicted=pred_class, Actual=df$group))

# =========================
# 12️⃣ Save Risk Scores
# =========================
df$risk_score <- pred

write.csv(df,
          "~/work/r_project_ishu/validation_blood/dataset-2/results/GSE15932_RiskScores_clean_april2.csv",
          row.names=FALSE)

# =========================
# 13️⃣ Save ROC Curve Coordinates
# =========================
roc_df <- data.frame(
  specificity = roc_panel$specificities,
  sensitivity = roc_panel$sensitivities
)

write.csv(roc_df,
          "~/work/r_project_ishu/validation_blood/dataset-2/results/GSE15932_ROC_curve_clean_april2.csv",
          row.names=FALSE)

cat("\n✅ External validation completed successfully.\n")
