# ============================================================
# Venn Diagram: WGCNA Hub Genes vs DEGs
# ============================================================

library(tidyverse)
library(VennDiagram)
library(grid)

# -----------------------
# Load data
# -----------------------
hub <- read.csv("~/work/r_project_ishu/wgcna/results/feb11 results/WGCNA_hub_genes-feb11.csv")
deg <- read.csv("~/work/r_project_ishu/deg/deg_results/DEG_tumor_vs_adjacent_results.csv")

# -----------------------
# Extract gene lists
# -----------------------
hub_genes <- unique(hub$Gene)

deg_genes <- deg %>%
  filter(status != "Not Significant") %>%   # keep only significant DEGs
  pull(Gene) %>% unique()

# -----------------------
# Intersection
# -----------------------
hub_deg_intersect <- intersect(hub_genes, deg_genes)
write.table(hub_deg_intersect, "~/work/r_project_ishu/wgcna/results/feb11 results/Hub_DEG_overlap_genes-feb11.txt", quote=FALSE, row.names=FALSE)

cat("Hub genes:", length(hub_genes), "\n")
cat("DEGs:", length(deg_genes), "\n")
cat("Overlap genes:", length(hub_deg_intersect), "\n")

# -----------------------
# Venn diagram
# -----------------------
venn.plot <- venn.diagram(
  x = list(
    HubGenes = hub_genes,
    DEGs = deg_genes
  ),
  filename = NULL,
  fill = c("red", "green"),
  alpha = 0.5,
  cex = 2,
  
  # ❌ REMOVE default labels completely
  cat.cex = 0,
  
  margin = 0.1,
  main = "Overlap of WGCNA Hub Genes and DEGs"
)

grid.newpage()
grid.draw(venn.plot)

png("~/work/r_project_ishu/venn/deg_wgcna_venn/deg_wgcna_plots/WGCNA_DEG_Venn-APRIL4.png",
    width = 1200, height = 1000, res = 150)

grid.newpage()

# Draw Venn
grid.draw(venn.plot)

# Add labels INSIDE device
grid.text("DEGs",
          x = 0.28, y = 0.82,
          gp = gpar(fontsize = 16, fontface = "bold"))

grid.text("WGCNA Hub Genes",
          x = 0.72, y = 0.82,
          gp = gpar(fontsize = 16, fontface = "bold"))

dev.off()