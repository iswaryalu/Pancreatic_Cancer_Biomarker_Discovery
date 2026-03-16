# ======================================================
# GSE62452 BIOMARKER VALIDATION (5-GENE PANEL)
# ======================================================

rm(list=ls())
set.seed(123)

# ------------------------
# 1️⃣ Load libraries
# ------------------------
library(pROC)
library(ggpubr)
library(dplyr)

# ------------------------
# 2️⃣ Define final biomarkers
# ------------------------
final_genes <- c("AHNAK2","TRIM29","SLC6A14","ITGB4","CTSE")

# ------------------------
# 3️⃣ Load expression matrix and metadata
# ------------------------
expr <- read.csv("~/work/r_project_ishu/validation_tumor/data/GSE62452_RAW_normalized_expression_matrix_pancreas-tumor-val.csv", 
                 row.names = 1, check.names = FALSE)

meta <- read.csv("~/work/r_project_ishu/validation_tumor/data/GSE62452_RAW_metadata_pancreas_tissue_only--tumor-val.csv", 
                 row.names = 1, check.names = FALSE)

# ------------------------
# 4️⃣ Extract only the samples present in metadata
# ------------------------
expr <- expr[, rownames(meta)]  # columns = sample IDs

# ------------------------
# 5️⃣ Recode tissue as group factor (Tumor vs Control)
# ------------------------
meta$group <- factor(ifelse(meta$tissue == "tumor", "Tumor", "Control"),
                     levels = c("Control", "Tumor"))

# Sanity check
table(meta$group)

# ------------------------
# 6️⃣ Extract biomarker genes present in matrix
# ------------------------
common_genes <- intersect(final_genes, rownames(expr))
expr_biomarker <- t(expr[common_genes, ])  # samples as rows
expr_biomarker <- as.data.frame(expr_biomarker)

# Combine group and expression
df <- data.frame(group = meta$group, expr_biomarker, check.names = FALSE)
cat("Genes used for validation:", common_genes, "\n")

# ------------------------
# 7️⃣ Boxplots of expression
# ------------------------
for(g in common_genes){
  p <- ggboxplot(df, x="group", y=g, color="group", add="jitter") +
    stat_compare_means(method="wilcox.test") +
    ggtitle(paste("GSE62452 Tumor Expression:", g))
  print(p)
}

# ------------------------
# 8️⃣ Single-gene ROC curves
# ------------------------
roc_table <- data.frame(Gene=character(), AUC=numeric())

for(g in common_genes){
  roc_obj <- roc(df$group, df[,g])
  plot(roc_obj, main=paste("ROC:", g))
  roc_table <- rbind(roc_table, data.frame(Gene=g, AUC=auc(roc_obj)))
}

write.csv(roc_table, "~/work/r_project_ishu/validation_tumor/results/GSE62452_SingleGene_AUC-feb18.csv", row.names=FALSE)
# Print single-gene AUCs
cat("Single-gene AUCs for GSE62452:\n")
print(roc_table)
# ------------------------
# 9️⃣ Multi-gene logistic regression ROC
# ------------------------
glm_model <- glm(group ~ ., data=df, family="binomial")
pred <- predict(glm_model, type="response")

roc_panel <- roc(df$group, pred)
# Open PNG device
png("~/work/r_project_ishu/validation_tumor/plots/GSE62452_Multigene_AUC_plot-feb18.png",
    width = 6, height = 4, units = "in", res = 300)

# Plot ROC
plot(roc_panel,
     col="blue",
     lwd=3,
     main="GSE62452 5-Gene Logistic ROC",
     print.auc=TRUE,
     print.auc.cex=1.2,
     print.auc.y=0.2)

# Close device
dev.off()

cat("Multi-gene AUC:", auc(roc_panel), "\n")
# Save ROC curve values
roc_df <- data.frame(
  specificity = roc_panel$specificities,
  sensitivity = roc_panel$sensitivities
)
write.csv(roc_df, "~/work/r_project_ishu/validation_tumor/results/GSE62452_Logistic_ROC_curve-feb18.csv", row.names=FALSE)

# ------------------------
# 🔟 Risk score output
# ------------------------
df$risk_score <- pred
write.csv(df, "~/work/r_project_ishu/validation_tumor/results/GSE62452_Biomarker_RiskScore-feb18.csv", row.names=FALSE)

# ------------------------
# 11️⃣ Confidence intervals for AUC
# ------------------------
ci_auc <- ci.auc(roc_panel)
ci_auc_boot <- ci.auc(roc_panel, method="bootstrap", boot.n=2000)
print(ci_auc)
print(ci_auc_boot)

# ------------------------
# 12️⃣ Best threshold metrics
# ------------------------
coords_best <- coords(roc_panel, "best", ret=c("threshold","sensitivity","specificity","ppv","npv"))
print(coords_best)
library(ggplot2)

# Reorder genes by AUC for plotting (highest first)
roc_table$Gene <- factor(roc_table$Gene, levels = roc_table$Gene[order(roc_table$AUC)])

# Create horizontal bar plot
p <- ggplot(roc_table, aes(x = Gene, y = AUC, fill = Gene)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(AUC, 3)), hjust = -0.1, size = 4) +
  coord_flip() +  # horizontal bars
  ylim(0, 1) +
  labs(title = "Single-Gene AUCs (GSE62452)", x = "Gene", y = "AUC") +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 12)) +
  scale_fill_brewer(palette = "Set2")

# Display the plot
print(p)

# Save the plot as PNG
ggsave("~/work/r_project_ishu/validation_tumor/plots/GSE62452_SingleGene_AUC_plot-feb18.png", plot = p, width = 6, height = 4, dpi = 300)

