# ============================================================
# Integration of DEGs, Hub Genes, and Significant Module Genes
# ============================================================

library(VennDiagram)
library(grid)       # for grid.draw
library(tidyverse)

# -----------------------
# 1. Load data
# -----------------------
deg <- read.csv("~/work/r_project_ishu/plots/DEG_tumor_vs_adjacent_results.csv")   # DEG table
sigModuleGenes <- read.csv("~/work/WGCNA_significant_module_genes_full_stats.csv")
hub_genes <- read.csv("~/work/WGCNA_hub_genes_full_stats.csv")

# -----------------------
# 2. Filter genes
# -----------------------
deg_filtered <- deg %>% filter(status != "Not Significant") %>% rename(Gene = 1)
sigModuleGenes_filtered <- sigModuleGenes %>% filter(SignificantModule == TRUE)
hub_genes_filtered <- hub_genes %>% filter(HubGene == TRUE & SignificantModule == TRUE)

# -----------------------
# 3. Find intersections with stats
# -----------------------
get_intersection_stats <- function(intersect_genes, df) {
  df %>% filter(Gene %in% intersect_genes)
}

# Hub ∩ DEGs
hub_deg_intersect_genes <- intersect(hub_genes_filtered$Gene, deg_filtered$Gene)
hub_deg_intersect <- get_intersection_stats(hub_deg_intersect_genes, hub_genes_filtered)
write.csv(hub_deg_intersect, "hub_deg_intersect_full_stats.csv", row.names = FALSE)

# Module ∩ DEGs
module_deg_intersect_genes <- intersect(sigModuleGenes_filtered$Gene, deg_filtered$Gene)
module_deg_intersect <- get_intersection_stats(module_deg_intersect_genes, sigModuleGenes_filtered)
write.csv(module_deg_intersect, "module_deg_intersect_full_stats.csv", row.names = FALSE)

# Hub ∩ Module ∩ DEGs
hub_module_deg_intersect_genes <- Reduce(intersect, list(
  hub_genes_filtered$Gene,
  sigModuleGenes_filtered$Gene,
  deg_filtered$Gene
))
hub_module_deg_intersect <- get_intersection_stats(hub_module_deg_intersect_genes, hub_genes_filtered)
write.csv(hub_module_deg_intersect, "hub_module_deg_intersect_full_stats.csv", row.names = FALSE)

cat("Number of hub genes overlapping DEGs:", nrow(hub_deg_intersect), "\n")
cat("Number of module genes overlapping DEGs:", nrow(module_deg_intersect), "\n")
cat("Number of hub genes in module overlapping DEGs:", nrow(hub_module_deg_intersect), "\n")

# -----------------------
# 4. Venn diagrams
# -----------------------
# Full Venn
venn_full_grob <- venn.diagram(
  x = list(HubGenes = hub_genes_filtered$Gene,
           ModuleGenes = sigModuleGenes_filtered$Gene,
           DEGs = deg_filtered$Gene),
  filename = NULL,   # NULL returns grob for drawing
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Integration of Hub Genes, Module Genes, and DEGs"
)
grid.newpage()
grid.draw(venn_full_grob)
# Save to PNG
venn.diagram(
  x = list(HubGenes = hub_genes_filtered$Gene,
           ModuleGenes = sigModuleGenes_filtered$Gene,
           DEGs = deg_filtered$Gene),
  filename = "WGCNA_DEG_Venn.png",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Integration of Hub Genes, Module Genes, and DEGs"
)

# DEG ∩ Hub genes
venn_hub_deg_grob <- venn.diagram(
  x = list(DEGs = deg_filtered$Gene,
           HubGenes = hub_genes_filtered$Gene),
  filename = NULL,
  fill = c("green", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "DEGs ∩ Hub Genes"
)
grid.newpage()
grid.draw(venn_hub_deg_grob)
png("DEG_HubGenes_Venn.png", width = 1000, height = 800, res = 150)
grid.draw(venn_hub_deg_grob)
dev.off()

# DEG ∩ Significant Module genes
venn_module_deg_grob <- venn.diagram(
  x = list(DEGs = deg_filtered$Gene,
           ModuleGenes = sigModuleGenes_filtered$Gene),
  filename = NULL,
  fill = c("green", "blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "DEGs ∩ Module Genes"
)
grid.newpage()
grid.draw(venn_module_deg_grob)
png("DEG_ModuleGenes_Venn.png", width = 1000, height = 800, res = 150)
grid.draw(venn_module_deg_grob)
dev.off()
