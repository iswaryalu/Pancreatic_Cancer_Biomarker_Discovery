```r
############################################################
# MASTER WGCNA PIPELINE (Tumor vs Adjacent Non-Tumor)
# Author: Iswarya (Computational Biology)
# Purpose: Paper-ready WGCNA with hub gene extraction
############################################################

# ================================
# 0. SETUP
# ================================
setwd("~/work/WGCNA_final")

pkgs <- c("WGCNA","tidyverse","matrixStats","ggplot2","pheatmap")
for(p in pkgs){ if(!require(p, character.only=TRUE)) install.packages(p, dependencies=TRUE) }

library(WGCNA)
library(tidyverse)
library(matrixStats)
library(ggplot2)
library(pheatmap)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ================================
# 1. LOAD DATA
# ================================
expr <- read.csv("~/work/r_project_ishu/scripts/normalized_expression_matrix.csv", row.names = 1)
meta <- read.csv("~/work/r_project_ishu/scripts/meta_clean.csv", row.names = 1)

# Clean sample names if needed
colnames(expr) <- gsub("_.*", "", colnames(expr))

# ================================
# 2. VARIANCE FILTERING (IMPORTANT)
# ================================
expr_mat <- as.matrix(expr)
gene_var <- rowVars(expr_mat)
expr_mat <- expr_mat[gene_var > quantile(gene_var, 0.25), ]
expr <- as.data.frame(expr_mat)
cat("Genes after filtering:", nrow(expr), "\n")

# ================================
# 3. SUBSET TUMOR VS ADJACENT
# ================================
meta_sub <- meta %>% filter(tissue %in% c("Tumor","adjacent_non_tumor"))

# Transpose expression (samples x genes)
datExpr <- t(expr[, rownames(meta_sub)])
datExpr <- as.data.frame(datExpr)
datExpr[] <- lapply(datExpr, as.numeric)

# Align metadata
rownames(meta_sub) <- rownames(datExpr)
stopifnot(all(rownames(meta_sub)==rownames(datExpr)))

# Trait matrix
traits_num <- data.frame(Tumor = ifelse(meta_sub$tissue=="Tumor",1,0))
rownames(traits_num) <- rownames(meta_sub)

# ================================
# 4. QUALITY CONTROL
# ================================
gsg <- goodSamplesGenes(datExpr, verbose=3)
if(!gsg$allOK){
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Sample clustering
sampleTree <- hclust(dist(datExpr), method="average")
png("01_sample_clustering.png",1200,800)
plot(sampleTree, main="Sample clustering")
abline(h=150,col="red")
dev.off()

# PCA
pca <- prcomp(datExpr, scale.=TRUE)
pca_df <- data.frame(pca$x[,1:2], Tissue=meta_sub$tissue)
png("02_PCA.png",1200,800)
ggplot(pca_df, aes(PC1,PC2,color=Tissue))+geom_point(size=3)+theme_minimal()
dev.off()

# ================================
# 5. SOFT THRESHOLD POWER
# ================================
power <- c(1:10, seq(12,50,2))
sft <- pickSoftThreshold(datExpr, powerVector=power, networkType="signed")

png("03_soft_threshold.png",1200,800)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Power", ylab="Scale Free R^2", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=power, col="red")
abline(h=0.8,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Power", ylab="Mean Connectivity", type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, col="red")
dev.off()

soft_power <- 14  # Change if needed

# ================================
# 6. NETWORK CONSTRUCTION
# ================================
bwnet <- blockwiseModules(datExpr,
                          power=soft_power,
                          TOMType="signed",
                          mergeCutHeight=0.25,
                          numericLabels=FALSE,
                          randomSeed=1234,
                          verbose=3)

moduleColors <- bwnet$colors
MEs <- bwnet$MEs

# Dendrogram
png("04_module_dendrogram.png",1400,900)
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  moduleColors[bwnet$blockGenes[[1]]],
  "Modules",
  dendroLabels=FALSE,
  hang=0.03
)
dev.off()

# ================================
# 7. MODULE-TRAIT HEATMAP
# ================================
moduleTraitCor <- cor(MEs, traits_num)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

textMatrix <- paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue,1),")",sep="")

png("05_module_trait_heatmap.png",1600,1200)
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels=colnames(traits_num),
               yLabels=colnames(MEs),
               ySymbols=colnames(MEs),
               colorLabels=FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix)
dev.off()

# ================================
# 8. MM & GS CALCULATION
# ================================
all_genes <- colnames(datExpr)
allModulesTable <- data.frame()

for(module in colnames(MEs)){
  module_color <- gsub("ME","",module)
  inModule <- moduleColors == module_color
  modGenes <- all_genes[inModule]
  
  MM <- cor(datExpr[,inModule], MEs[,module])
  MM_pval <- corPvalueStudent(MM, nrow(datExpr))
  GS <- cor(datExpr[,inModule], traits_num$Tumor)
  GS_pval <- corPvalueStudent(GS, nrow(datExpr))
  
  df <- data.frame(Gene=modGenes, Module=module_color,
                   MM=as.vector(MM), MM_pval=as.vector(MM_pval),
                   GS_Tumor=as.vector(GS), GS_pval=as.vector(GS_pval))
  allModulesTable <- rbind(allModulesTable, df)
}

# 9. Identify Significant Modules
# -----------------------
sigModules <- rownames(moduleTraitCor)[
  abs(moduleTraitCor[, "Tumor"]) > 0.5 &
    moduleTraitPvalue[, "Tumor"] < 0.05
]

cat("Significant modules:", sigModules, "\n")

# -----------------------
# 10. Identify Hub Genes (Paper Standard)
# -----------------------
allModulesTable$SignificantModule <- allModulesTable$Module %in% gsub("ME","",sigModules)

hubGenesTable <- allModulesTable %>%
  filter(SignificantModule == TRUE,
         abs(MM) > 0.8,
         abs(GS_Tumor) > 0.2)

# ================================
# 11. SAVE OUTPUTS
# ================================
write.csv(allModulesTable, "WGCNA_all_genes_MM_GS-feb11.csv", row.names=FALSE)
write.csv(hubGenesTable, "WGCNA_hub_genes-feb11.csv", row.names=FALSE)
write.table(sigModules, "WGCNA_significant_modules-feb11.txt", quote=FALSE, row.names=FALSE)
write.csv(allModulesTable %>% select(Gene,Module), "WGCNA_gene_module_mapping-feb11.csv", row.names=FALSE)

cat("\nWGCNA COMPLETE: Paper-ready outputs generated\n")
```
