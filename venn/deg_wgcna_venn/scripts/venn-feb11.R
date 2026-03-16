# ============================================================
# Venn Diagram: WGCNA Hub Genes vs DEGs
# ============================================================

library(tidyverse)
library(VennDiagram)
library(grid)

# -----------------------
# Load data
# -----------------------
hub <- read.csv("WGCNA_hub_genes-feb11.csv")
deg <- read.csv("~/work/r_project_ishu/plots/DEG_tumor_vs_adjacent_results.csv")

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
write.table(hub_deg_intersect, "Hub_DEG_overlap_genes-feb11.txt", quote=FALSE, row.names=FALSE)

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
  cat.cex = 1.5,
  main = "Overlap of WGCNA Hub Genes and DEGs"
)

grid.newpage()
grid.draw(venn.plot)

# Save figure
png("WGCNA_DEG_Venn-feb11.png", width=1200, height=1000, res=150)
grid.draw(venn.plot)
dev.off()
