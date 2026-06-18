# ============================================================
# WGCNA for Tumor vs Adjacent Non-Tumor
# ============================================================

# -----------------------
# 0. Setup
# -----------------------
setwd("~/work")

# Install & load packages
if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("factoextra", quietly = TRUE)) install.packages("factoextra")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

library(WGCNA)
library(tidyverse)
library(matrixStats)
library(devtools)
library(factoextra)
library(gridExtra)

# Install CorLevelPlot (optional)
devtools::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# -----------------------
# 1. Load data
# -----------------------
expr <- read.csv("~/work/r_project_ishu/scripts/normalized_expression_matrix.csv", row.names = 1)
meta <- read.csv("~/work/r_project_ishu/scripts/meta_clean.csv", row.names = 1)

# Clean sample names
colnames(expr) <- gsub("_.*", "", colnames(expr))

# -----------------------
# 2. Subset for Tumor vs Adjacent Non-Tumor
# -----------------------
meta_sub <- meta %>% filter(tissue %in% c("Tumor", "adjacent_non_tumor"))
datExpr_sub <- t(expr[, rownames(meta_sub)])
datExpr_sub <- as.data.frame(datExpr_sub)
datExpr_sub[] <- lapply(datExpr_sub, as.numeric)

# Ensure rownames match
rownames(meta_sub) <- rownames(datExpr_sub)
stopifnot(all(rownames(meta_sub) == rownames(datExpr_sub)))

# -----------------------
# 3. Check good samples and genes
# -----------------------
gsg <- goodSamplesGenes(datExpr_sub, verbose = 3)
if (!gsg$allOK) { stop("Some genes or samples are problematic. Inspect gsg object.") }
gsg$allOK

# -----------------------
# 4. Sample clustering & PCA
# -----------------------
# Sample clustering
sampleTree <- hclust(dist(datExpr_sub), method = "average")
png("sample_tree.png", width = 1000, height = 600)
plot(sampleTree, main = "Sample clustering (Tumor vs Adjacent Non-Tumor)", sub = "", xlab = "", cex = 0.6)
abline(h = 150, col = "red")
dev.off()

# PCA
pca <- prcomp(datExpr_sub, scale. = TRUE)
pca_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Tissue = meta_sub$tissue)

png("PCA_plot.png", width = 1000, height = 800)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Tumor vs Adjacent Non-Tumor",
       x = paste0("PC1 (", pca_var[1], "% variance)"),
       y = paste0("PC2 (", pca_var[2], "% variance)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# -----------------------
# 5. Soft-thresholding
# -----------------------
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
sft <- pickSoftThreshold(datExpr_sub, powerVector = power, networkType = "signed", verbose = 5)
sft.data <- sft$fitIndices

# Plot scale-free topology
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() + geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Soft threshold (Power)", y = "Scale-free topology model fit, signed R^2") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  labs(x = "Soft threshold (Power)", y = "Mean Connectivity") +
  theme_classic()
combined_plot <- arrangeGrob(a1, a2, nrow = 2)
ggsave("soft_thresholding_plot_tumor_vs_adj.pdf", combined_plot, width = 8, height = 10)
# Save the first plot (Scale-free topology)
ggsave(
  "soft_thresholding_scale_free_tumor_vs_adj.png",
  plot = a1,
  width = 6, height = 5, dpi = 300
)

# Save the second plot (Mean connectivity)
ggsave(
  "soft_thresholding_mean_connectivity_tumor_vs_adj.png",
  plot = a2,
  width = 7, height = 5, dpi = 300
)

# Choose soft-threshold
soft_power <- 14  # adjust based on your plot

# -----------------------
# 6. Network construction & module detection
# -----------------------
bwnet <- blockwiseModules(datExpr_sub,
                          maxBlockSize = 20000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

# Module dendrogram
png("module_dendrogram_tumor_vs_adj.png", width = 1200, height = 800, res = 150)
plotDendroAndColors(bwnet$dendrograms[[1]],
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

# -----------------------
png("module_trait_heatmap_tumor_vs_adj.png", width = 2000, height = 1600, res = 200)

par(mar = c(8, 10, 4, 6))  # increase bottom, left, top, right margins

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traits_num),
               yLabels = colnames(bwnet$MEs),
               ySymbols = colnames(bwnet$MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               cex.lab = 1.2,
               zlim = c(-1, 1),
               main = "Module–Tissue Relationships")

dev.off()

# -----------------------
# 8. Extract module/hub genes with full stats
# -----------------------
all_genes <- colnames(datExpr_sub)
all_modules <- bwnet$colors
module_names <- colnames(bwnet$MEs)
traits <- colnames(traits_num)

trait_of_interest <- "Tumor"  # or whichever column in traits_num
cor_threshold <- 0.5
pval_threshold <- 0.05

sigModules <- module_names[
  (abs(moduleTraitCor[, trait_of_interest]) > cor_threshold) & 
    (moduleTraitPvalue[, trait_of_interest] < pval_threshold)
]


# Initialize tables
allModulesTable <- data.frame()
sigModuleGenesTable <- data.frame()
hubGenesTable <- data.frame()

for (module in module_names) {
  module_color <- gsub("ME", "", module)
  inModule <- (all_modules == module_color)
  modGenes <- all_genes[inModule]
  
  MM <- cor(datExpr_sub[, inModule], bwnet$MEs[, module], use = "p")
  MM_pval <- corPvalueStudent(MM, nrow(datExpr_sub))
  
  GS <- sapply(traits, function(trait) cor(datExpr_sub[, inModule], traits_num[, trait], use = "p"))
  GS_pval <- sapply(traits, function(trait) corPvalueStudent(GS[, trait], nrow(datExpr_sub)))
  
  df <- data.frame(
    Gene = modGenes,
    Module = module_color,
    MM = as.vector(MM),
    MM_pval = as.vector(MM_pval),
    stringsAsFactors = FALSE
  )
  
  for (trait in traits) {
    df[[paste0("GS_", trait)]] <- GS[, trait]
    df[[paste0("pval_GS_", trait)]] <- GS_pval[, trait]
  }
  
  df$SignificantModule <- ifelse(module %in% sigModules, TRUE, FALSE)
  df$HubGene <- (abs(df$MM) > 0.8 & apply(abs(GS) > 0.2, 1, any))
  
  allModulesTable <- rbind(allModulesTable, df)
  sigModuleGenesTable <- rbind(sigModuleGenesTable, df[df$SignificantModule,])
  hubGenesTable <- rbind(hubGenesTable, df[df$HubGene,])
}

# -----------------------
# 9. Save CSVs
# -----------------------
gene_module <- allModulesTable %>% select(Gene, Module)
write.csv(gene_module, "WGCNA_gene_module_mapping_1.csv", row.names = FALSE)




write.csv(moduleTraitCor, "module_trait_correlation_tumor_vs_adj.csv", row.names = TRUE)
write.csv(moduleTraitPvalue, "module_trait_pvalues_tumor_vs_adj.csv", row.names = TRUE)

write.csv(allModulesTable, "WGCNA_all_modules_full_stats.csv", row.names = FALSE)
write.csv(sigModuleGenesTable, "WGCNA_significant_module_genes_full_stats.csv", row.names = FALSE)
write.csv(hubGenesTable, "WGCNA_hub_genes_full_stats.csv", row.names = FALSE)

cat("\n✅ All CSVs saved successfully.\n")
