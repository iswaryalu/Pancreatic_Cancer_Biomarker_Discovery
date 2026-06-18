# -----------------------------------------------------------
# 1. Data Correction, Normalization, and Annotation
# -----------------------------------------------------------
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GSEABase")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library(limma)
library(clusterProfiler)
library(GSEABase)
library(tidyverse)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
# Install ggrepel if not installed yet
if(!require(ggrepel)) install.packages("ggrepel")

# Load the library
library(ggrepel)


expr <- read.csv("~/work/r_project_ishu/scripts/normalized_expression_matrix.csv", row.names = 1)
meta <-read.csv("~/work/r_project_ishu/scripts/meta_clean.csv", row.names = 1)
dim(expr)
dim(meta)
colnames(expr) <- gsub("_.*", "", colnames(expr))
head(colnames(expr))
all(colnames(expr) %in% rownames(meta))
all(rownames(meta) %in% colnames(expr))


# -----------------------------------------------------------
# 2. DEG analysis using limma
# -----------------------------------------------------------

meta_dge <- meta %>% filter(tissue %in% c("Tumor", "adjacent_non_tumor"))
expr_dge <- expr[, rownames(meta_dge)]
stopifnot(all(colnames(expr_dge) == rownames(meta_dge)))

group <- factor(meta_dge$tissue, levels = c("adjacent_non_tumor", "Tumor"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_dge, design)
contrast.matrix <- makeContrasts(Tumor - adjacent_non_tumor, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg_results <- topTable(fit2, number = Inf, adjust.method = "BH")

# Apply cutoff: |log2FC| > 1, p < 0.05
deg_results <- deg_results %>%
  rownames_to_column("Gene") %>%
  mutate(status = case_when(
    P.Value < 0.05 & logFC > 1 ~ "Upregulated",
    P.Value < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

write.csv(deg_results, "DEG_tumor_vs_adjacent_results.csv", row.names = FALSE)
# -----------------------------------------------------------
# 3. Visualization with ggplot2
# -----------------------------------------------------------


# Volcano plot with gray points for non-significant genes
library(ggplot2)
ggplot(deg_results,
       aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "Volcano Plot: Tumor vs Adjacent Non-Tumor",
       x = "log2 Fold Change",
       y = "-log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5))

library(ggplot2)
library(ggrepel)

ggplot(deg_results, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_text_repel(
    data = subset(deg_results, P.Value < 0.01 & abs(logFC) > 2),  # label top DEGs only
    aes(label = Gene),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 15
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Tumor vs Adjacent Non-Tumor",
    x = "log2 Fold Change",
    y = "-log10 P-value"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# -----------------------------------------------------------
# 4. Functional Enrichment using GSEABase
# -----------------------------------------------------------

# Load packages
library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)

# Step 1: prepare ranked gene list
geneList <- fit2$coefficients[, 1]
names(geneList) <- rownames(fit2$coefficients)
geneList <- sort(geneList, decreasing = TRUE)
names(geneList) <- toupper(names(geneList))  # ensure uppercase

# Step 2: load hallmark gene sets
gene.sets <- getGmt("~/work/r_project_ishu/scripts/h.all.v2025.1.Hs.symbols.gmt")

# Step 3: make TERM2GENE correctly
term2gene <- stack(geneIds(gene.sets))
colnames(term2gene) <- c("gene", "term")      # stack gives gene first
term2gene <- term2gene[, c("term", "gene")]   # reorder columns to term,gene
term2gene$gene <- toupper(term2gene$gene)

# Step 4: check overlap
length(intersect(names(geneList), term2gene$gene))

# Step 5: run GSEA
gsea_res <- GSEA(
  geneList     = geneList,
  TERM2GENE    = term2gene,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
head(gsea_res@result[, c("ID", "Description", "NES", "p.adjust")])


library(enrichplot)
dotplot(gsea_res, showCategory = 15, split = ".sign") + 
  facet_grid(.~.sign)

gseaplot2(gsea_res, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
gseaplot2(gsea_res, geneSetID = "HALLMARK_KRAS_SIGNALING_UP")
library(ggplot2)
library(enrichplot)

# Dotplot with wider color range and clearer separation
dotplot(gsea_res, showCategory = 20, split = ".sign") + 
  facet_grid(. ~ .sign) +
  scale_color_gradientn(
    colors = c("blue", "skyblue", "pink", "red"),
    limits = c(0, 0.05),
    name = "Adjusted p-value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

min_p <- min(gsea_res@result$p.adjust, na.rm = TRUE)
max_p <- max(gsea_res@result$p.adjust, na.rm = TRUE)
head(min_p)
head(max_p)
ridgeplot(gsea_res, showCategory = 25, fill = "p.adjust") +
  scale_fill_gradientn(
    colors = c("red", "white", "blue"),
    name = "Adjusted p-value"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

range(gsea_res@result$p.adjust)
ridgeplot(gsea_res, showCategory = 25, fill = "p.adjust") +
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),
    trans = "log10",  # spread small p-values visually
    name = "Adjusted p-value (log10)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )
ridgeplot(gsea_res, showCategory = 25, fill = "NES") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = "NES"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )
ridge_p <- ridgeplot(gsea_res, showCategory = 20, fill = "NES") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, 
    name = "NES"
  ) +
  labs(title = "GSEA Ridge Plot (Top Pathways)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

dot_p <- dotplot(gsea_res, showCategory = 20) +
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
    trans = "log10",
    name = "Adjusted p-value"
  ) +
  labs(title = "GSEA Dot Plot (Top Pathways)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
combined_plot <- ridge_p + dot_p + plot_layout(widths = c(1.2, 1))
combined_plot
library(plot_layout)
install.packages("patchwork")
library(patchwork)
# Example: DEG table (adjust this to your object name)
# deg_results <- your_DEG_results_dataframe
colnames(deg_results)
top30_genes <- deg_results %>%
  arrange(adj.P.Val) %>%     # smallest adjusted p-values
  slice(1:30) %>%
  pull(Gene)                 # column name of gene symbols
expr_top30 <- expr_dge[top30_genes, ]
rownames(expr_dge) <- deg_results$Gene
expr_top30_scaled <- t(scale(t(expr_top30)))
library(pheatmap)

# Example if you have a sample group info data frame:
# annotation_col <- data.frame(Group = meta$tissue)
# rownames(annotation_col) <- colnames(expr_top30_scaled)

pheatmap(
  expr_top30_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "Top 30 DEGs Expression Heatmap"
)
colnames(expr_dge)
# Step 4: Prepare annotation based on tissue type
annotation_col <- data.frame(
  Tissue = meta_dge$tissue
)
rownames(annotation_col) <- rownames(meta_dge)

# Step 5: Define colors for annotation
ann_colors <- list(
  Tissue = c(Tumor = "firebrick3", adjacent_non_tumor = "steelblue3")
)

# Step 6: Plot heatmap
pheatmap(
  expr_top30_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize_row = 8,
  main = "Top 30 DEGs: Tumor (Red) vs Adjacent Non-Tumor (Blue)"
)





















# Basic volcano plot
volcano_plot <- ggplot(deg_results,
                       aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(title = "Volcano Plot: Tumor vs Adjacent Non-Tumor",
       x = "log2 Fold Change",
       y = "-log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("~/work/r_project_ishu/plots/Volcano_Plot_basic.png", volcano_plot, width = 7, height = 5, dpi = 300)

# Volcano plot with top DEGs labeled
volcano_label_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_text_repel(
    data = subset(deg_results, P.Value < 0.01 & abs(logFC) > 2),
    aes(label = Gene),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 50
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Tumor vs Adjacent Non-Tumor",
    x = "log2 Fold Change",
    y = "-log10 P-value"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("~/work/r_project_ishu/plots/Volcano_Plot_labeled.png", volcano_label_plot, width = 10, height = 8, dpi = 300)
library(pheatmap)

# Top30 scaled expression matrix
expr_top30_scaled <- t(scale(t(expr_top30)))

# Heatmap with annotation
annotation_col <- data.frame(
  Tissue = meta_dge$tissue
)
rownames(annotation_col) <- rownames(meta_dge)
ann_colors <- list(
  Tissue = c(Tumor = "firebrick3", adjacent_non_tumor = "steelblue3")
)

# Save heatmap as PNG
png("~/work/r_project_ishu/plots/Top30_DEGs_Heatmap.png", width = 1200, height = 1000, res = 150)
pheatmap(
  expr_top30_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize_row = 8,
  main = "Top 30 DEGs: Tumor vs Adjacent Non-Tumor"
)
dev.off()
# Dotplot
dot_p <- dotplot(gsea_res, showCategory = 20) +
  scale_color_gradientn(
    colors = c("blue", "white", "red"),
    trans = "log10",
    name = "Adjusted p-value"
  ) +
  labs(title = "GSEA Dot Plot (Top Pathways)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave("GSEA_dotplot.png", dot_p, width = 8, height = 6, dpi = 300)

# Ridge plot
ridge_p <- ridgeplot(gsea_res, showCategory = 20, fill = "NES") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, 
    name = "NES"
  ) +
  labs(title = "GSEA Ridge Plot (Top Pathways)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave("GSEA_ridgeplot.png", ridge_p, width = 10, height = 6, dpi = 300)

# Combined plot using patchwork
library(patchwork)
combined_plot <- ridge_p + dot_p + plot_layout(widths = c(1.2, 1))
ggsave("GSEA_combined_plot.png", combined_plot, width = 15, height = 6, dpi = 300)
top30_genes <- deg_results %>% arrange(adj.P.Val) %>% slice(1:30)
write.csv(top30_genes, "DEG_top30.csv", row.names = FALSE)

head(deg_results)
# Gene      logFC      P.Value   adj.P.Val   status
# A1CF      1.5        0.002     0.01       Upregulated
# A2M      -1.2        0.03      0.04       Downregulated
