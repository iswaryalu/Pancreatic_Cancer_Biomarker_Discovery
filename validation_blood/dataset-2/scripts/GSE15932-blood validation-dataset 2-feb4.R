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
meta$group <- ifelse(meta$Group %in% c("PC+Diabetes","PC only"), "Tumor",
                     ifelse(meta$Group %in% c("Diabetes only","Healthy"), "Control", NA))

meta <- meta[!is.na(meta$group), ]
meta$group <- factor(meta$group, levels=c("Control","Tumor"))

cat("Group distribution:\n")
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

dir.create("~/work/r_project_ishu/validation_blood/dataset-2/results/single_gene_ROC",
           showWarnings = FALSE)

for (gene in common_genes) {
  
  roc_single <- roc(df$group, df[[gene]])
  
  auc_val <- auc(roc_single)
  
  single_auc <- rbind(single_auc,
                      data.frame(Gene=gene, AUC=as.numeric(auc_val)))
  
  # Save individual ROC
  png(paste0("~/work/r_project_ishu/validation_blood/dataset-2/plots/single_gene_ROC/",
             gene, "_ROC.png"),
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
          "~/work/r_project_ishu/validation_blood/dataset-2/results/GSE15932_single_gene_AUC.csv",
          row.names=FALSE)

# =========================
# 9️⃣ Plot ROC
# =========================
# =========================
# 9️⃣ Save Multi-Gene ROC (Base)
# =========================
png("~/work/r_project_ishu/validation_blood/dataset-2/plots/GSE15932_multi_gene_ROC.png",
    width=800, height=700)

plot(roc_panel,
     col="black",
     lwd=3,
     main="GSE15932 Multi-Gene Logistic ROC")

text(0.6, 0.2,
     labels=paste("AUC =", round(auc(roc_panel),3)),
     cex=1.2)

dev.off()

# Overlay ROC curves
# =========================
png("~/work/r_project_ishu/validation_blood/dataset-2/plots/GSE15932_ROC_overlay.png",
    width=800, height=700)

plot(roc_panel, col="black", lwd=3,
     main="Multi vs Single Gene ROC")

cols <- rainbow(length(common_genes))

i <- 1
for (gene in common_genes) {
  roc_single <- roc(df$group, df[[gene]])
  plot(roc_single, add=TRUE, col=cols[i], lwd=2)
  i <- i + 1
}

legend("bottomright",
       legend=c("Multi-gene", common_genes),
       col=c("black", cols),
       lwd=2)

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
threshold <- opt_cut["threshold"]

pred_class <- ifelse(pred >= threshold, "Tumor", "Control")
pred_class <- factor(pred_class, levels=c("Control","Tumor"))

cat("\nConfusion Matrix:\n")
print(table(Predicted=pred_class, Actual=df$group))

# =========================
# 12️⃣ Save Risk Scores
# =========================
df$risk_score <- pred

write.csv(df,
          "~/work/r_project_ishu/validation_blood/dataset-2/results/GSE15932_RiskScores_clean.csv",
          row.names=FALSE)

# =========================
# 13️⃣ Save ROC Curve Coordinates
# =========================
roc_df <- data.frame(
  specificity = roc_panel$specificities,
  sensitivity = roc_panel$sensitivities
)

write.csv(roc_df,
          "~/work/r_project_ishu/validation_blood/dataset-2/results/GSE15932_ROC_curve_clean.csv",
          row.names=FALSE)

cat("\n✅ External validation completed successfully.\n")
