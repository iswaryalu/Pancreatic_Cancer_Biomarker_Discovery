# ============================
# 1. INSTALL & LOAD PACKAGES
# ============================
if (!require(VennDiagram)) install.packages("VennDiagram")
if (!require(ggvenn)) install.packages("ggvenn")
library(VennDiagram)
library(ggvenn)

# ============================
# 2. LOAD GENE LISTS
# ============================
svm  <- read.csv("~/work/r_project_ishu/ml/svm-rfe/results/SVM_RFE_biomarker_genes-feb18.csv")
rf   <- read.csv("RF_genes_importance_gt1.csv")
lasso <- read.csv("~/work/r_project_ishu/ml/lasso/results/LASSO_selected_biomarkers-feb17.csv")

# ============================
# 3. EXTRACT GENE NAMES (COLUMN = x)
# ============================
svm_genes   <- unique(svm$x)
rf_genes    <- unique(rf$Gene)
lasso_genes <- unique(lasso$x)

# ============================
# 4. CREATE LIST FOR VENN
# ============================
gene_sets <- list(
  SVM_RFE = svm_genes,
  Random_Forest = rf_genes,
  LASSO = lasso_genes
)

# ============================
# 5. CLASSIC VENN DIAGRAM
# ============================
venn.diagram(
  x = gene_sets,
  filename = "Biomarker_VennDiagram.png",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0,
  main = "Overlap of Biomarker Genes from ML Methods"
)

# ============================
# 6. MODERN PUBLICATION STYLE VENN
# ============================
ggvenn(
  gene_sets,
  fill_color = c("red", "blue", "green"),
  stroke_size = 0.5,
  set_name_size = 6
)

ggsave("~/work/r_project_ishu/venn/ml_venn_plots/Biomarker_VennDiagram_ggvenn-feb18.png", width=6, height=6, dpi=300)
common_all <- Reduce(intersect, gene_sets)
print(common_all)

write.csv(common_all, "~/work/r_project_ishu/venn/biomarker_results/Common_Biomarker_Genes_All3-feb18.csv", row.names = FALSE)
