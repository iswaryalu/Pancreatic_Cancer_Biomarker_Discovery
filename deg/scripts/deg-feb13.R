# ==========================================
# FIGURE 2 – DEG + Heatmap + GSEA
# ==========================================

rm(list = ls())

library(limma)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(pheatmap)

# ------------------------------
# 1. Load Data
# ------------------------------

expr <- read.csv("~/work/r_project_ishu/normalized_expression_matrix.csv",
                 row.names = 1, check.names = FALSE)

meta <- read.csv("~/work/r_project_ishu/meta_clean.csv",
                 row.names = 1)

colnames(expr) <- gsub("_.*", "", colnames(expr))

meta_dge <- meta %>%
  filter(tissue %in% c("Tumor", "adjacent_non_tumor"))

expr_dge <- expr[, rownames(meta_dge)]

# ------------------------------
# 2. limma DEG
# ------------------------------

group <- factor(meta_dge$tissue,
                levels = c("adjacent_non_tumor", "Tumor"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_dge, design)

contrast.matrix <- makeContrasts(Tumor - adjacent_non_tumor,
                                 levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(fit2, number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Gene")

# DEG classification
deg_results <- deg_results %>%
  mutate(status = case_when(
    adj.P.Val < 0.05 & logFC > 1  ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

write.csv(deg_results,
          "~/work/r_project_ishu/deg/feb17/DEG_results_feb17.csv",
          row.names = FALSE)
# =========================================================
# PART 0 — Volcano Plot
# =========================================================

library(ggrepel)

volcano_plot <- ggplot(deg_results, 
                       aes(x = logFC, 
                           y = -log10(adj.P.Val), 
                           color = status)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(title = "Volcano Plot (Tumor vs Adjacent)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")

ggsave("~/work/r_project_ishu/deg/feb17/Volcano_plot_feb17.png",
       volcano_plot,
       width = 7, height = 6)

# =========================================================
# PART A — Heatmap of Top 30 DEGs
# =========================================================

top30 <- deg_results %>%
  arrange(adj.P.Val) %>%
  filter(status != "Not Significant") %>%
  head(30)

top30_genes <- top30$Gene

expr_top30 <- expr_dge[top30_genes, ]

# Z-score scaling
expr_scaled <- t(scale(t(expr_top30)))

annotation_col <- data.frame(
  Tissue = meta_dge$tissue
)
rownames(annotation_col) <- rownames(meta_dge)

pheatmap(expr_scaled,
         annotation_col = annotation_col,
         show_colnames = FALSE,
         fontsize_row = 8,
         clustering_method = "complete",
         filename = "~/work/r_project_ishu/deg/feb17/Top30_DEG_heatmap_feb17.png")

library(scales)

# Prepare ranked gene list (ALL genes, no cutoff)
gene_list <- deg_results$logFC
names(gene_list) <- toupper(deg_results$Gene)

# Convert SYMBOL to ENTREZ
gene_df_all <- bitr(names(gene_list),
                    fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Hs.eg.db)

gene_list <- gene_list[gene_df_all$SYMBOL]
names(gene_list) <- gene_df_all$ENTREZID

gene_list <- sort(gene_list, decreasing = TRUE)

# ------------------------------
# GSEA KEGG
# ------------------------------

gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

# Sort pathways by adjusted p-value
gsea_kegg@result <- gsea_kegg@result %>%
  arrange(p.adjust)

write.csv(as.data.frame(gsea_kegg),
          "~/work/r_project_ishu/deg/deg_results/GSEA_KEGG_march5.csv")

# ------------------------------
# Ridge Plot (Top 30 pathways)
# ------------------------------

ridge_plot <- ridgeplot(
  gsea_kegg,
  showCategory = 30,
  fill = "p.adjust"
) +
  scale_fill_gradientn(
    colours = c("red", "purple", "blue"),
    name = "p.adjust",
    trans = "log10"
  ) +
  labs(
    title = "GSEA KEGG Pathway Enrichment",
    x = "NES"
  ) +
  theme_bw(base_size = 12)

ridge_plot
ggsave("~/work/r_project_ishu/deg/deg_plots/GSEA_ridgeplot_march5.png",
       ridge_plot,
       width = 8,
       height = 6,
       dpi = 300)
# =========================================================
# PART C — Export Genes for CMap
# =========================================================

library(dplyr)

deg <- deg_results

# Remove NA genes
deg <- deg %>% filter(!is.na(Gene))

# UP genes
up_genes <- deg %>%
  filter(status == "Upregulated") %>%
  arrange(desc(logFC)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  head(100)

# DOWN genes
down_genes <- deg %>%
  filter(status == "Downregulated") %>%
  arrange(logFC) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  head(100)

# Save for CMap
write.table(up_genes$Gene,
            "~/work/r_project_ishu/cmap_up_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(down_genes$Gene,
            "~/work/r_project_ishu/cmap_down_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("CMap gene lists exported successfully!\n")
unique(up_genes$Gene)
cat("Up genes:", nrow(up_genes), "\n")
cat("Down genes:", nrow(down_genes), "\n")