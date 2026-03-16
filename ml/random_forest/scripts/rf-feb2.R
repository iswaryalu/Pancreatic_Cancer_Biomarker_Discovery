rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

# =========================
# 1. LOAD LIBRARIES
# =========================
packages <- c("randomForest","ggplot2","dplyr")
for(p in packages){
  if(!require(p, character.only=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

cat("Libraries loaded\n")

# =========================
# 2. LOAD DATA
# =========================
expr  <- read.csv("normalized_expression_matrix.csv", row.names = 1, check.names = FALSE)
genes <- read.csv("Hub_DEG_overlap_genes-feb11.txt")
meta  <- read.csv("meta_clean.csv", stringsAsFactors = FALSE)

gene_list <- unique(genes$Gene)

# =========================
# 3. PREPARE DATA
# =========================
expr_filt <- expr[rownames(expr) %in% gene_list, ]
expr_t <- as.data.frame(t(expr_filt))
expr_t$SampleID <- gsub("_.*", "", rownames(expr_t))

merged <- merge(meta, expr_t, by="SampleID")

# Tumor vs Adjacent only
merged <- merged[merged$tissue %in% c("Tumor","adjacent_non_tumor"), ]
merged$tissue <- factor(merged$tissue,
                        levels=c("adjacent_non_tumor","Tumor"))

ml_data <- merged[, c("tissue", gene_list)]

x <- data.matrix(ml_data[, gene_list])
y <- ml_data$tissue

# Remove NA
idx <- complete.cases(x)
x <- x[idx,]
y <- y[idx]

cat("Samples:", nrow(x), " Genes:", ncol(x), "\n")

# =========================
# 4. RANDOM FOREST MODEL
# =========================
rf_model <- randomForest(
  x = x,
  y = y,
  ntree = 2000,
  importance = TRUE
)

cat("OOB Error:",
    rf_model$err.rate[2000,"OOB"], "\n")

# =================================================
# 📊 PANEL C — ERROR RATE VS TREES (PAPER FIGURE)
# =================================================
png("RF_error_vs_trees.png", width=900, height=700)

plot(rf_model,
     main="Error Rate vs Number of Trees",
     lwd=2,
     col=c("black","red","green"))

legend("topright",
       legend=c("OOB","Adjacent","Tumor"),
       col=c("black","red","green"),
       lwd=2,
       bty="n")

dev.off()

# =================================================
# 📊 FEATURE IMPORTANCE (Gini)
# =================================================
imp <- importance(rf_model, type=2)

imp_df <- data.frame(
  Gene = rownames(imp),
  Importance = imp[,1]
)


# Sort by importance
imp_df <- imp_df[order(-imp_df$Importance), ]

# Save full ranking
write.csv(imp_df, "RF_gene_importance_ranking.csv", row.names=FALSE)

# =================================================
# ⭐ SELECT GENES WITH IMPORTANCE > 1 (PAPER METHOD)
# =================================================
rf_selected <- imp_df[imp_df$Importance > 1, ]

cat("Genes with importance > 1:", nrow(rf_selected), "\n")

# Save selected genes
write.csv(rf_selected,
          "RF_genes_importance_gt1.csv",
          row.names = FALSE)

# =================================================
# 📊 PANEL D — TOP 30 GINI IMPORTANCE
# =================================================
top30 <- imp_df[1:30, ]

p <- ggplot(top30,
            aes(x=Importance,
                y=reorder(Gene, Importance))) +
  geom_point(size=3) +
  theme_bw() +
  labs(title="Top 30 Genes by Random Forest",
       x="Mean Decrease in Gini",
       y="Feature")

ggsave("RF_top30_gini.png", p, width=6, height=8, dpi=300)

# =================================================
# 📄 EXPORT TOP GENES FOR INTERSECTION
# =================================================
top_rf_genes <- imp_df$Gene[1:30]
write.csv(top_rf_genes,
          "RF_top30_genes_for_intersection.csv",
          row.names=FALSE)

cat("\n✅ Random Forest analysis completed.\n")
