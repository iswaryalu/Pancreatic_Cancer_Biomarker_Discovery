# ==========================================
# 02. GSEA Analysis (Hallmark)
# ==========================================

library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)
library(ggplot2)

# ------------------------------
# 1. Load DEG fit object OR re-run ranking
# ------------------------------

expr <- read.csv("~/work/r_project_ishu/scripts/normalized_expression_matrix.csv",
                 row.names = 1, check.names = FALSE)

meta <- read.csv("~/work/r_project_ishu/scripts/meta_clean.csv",
                 row.names = 1)

colnames(expr) <- gsub("_.*", "", colnames(expr))

meta_dge <- meta %>%
  filter(tissue %in% c("Tumor", "adjacent_non_tumor"))

expr_dge <- expr[, rownames(meta_dge)]

group <- factor(meta_dge$tissue,
                levels = c("adjacent_non_tumor", "Tumor"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_dge, design)
contrast.matrix <- makeContrasts(Tumor - adjacent_non_tumor,
                                 levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# ------------------------------
# 2. Create Ranked Gene List
# ------------------------------

geneList <- fit2$t[, 1]   # moderated t-statistic
names(geneList) <- rownames(fit2)
geneList <- sort(geneList, decreasing = TRUE)
names(geneList) <- toupper(names(geneList))

# ------------------------------
# 3. Load Hallmark Gene Sets
# ------------------------------

gene.sets <- getGmt("~/work/r_project_ishu/scripts/h.all.v2025.1.Hs.symbols.gmt")

term2gene <- stack(geneIds(gene.sets))
colnames(term2gene) <- c("gene", "term")
term2gene <- term2gene[, c("term", "gene")]
term2gene$gene <- toupper(term2gene$gene)

# ------------------------------
# 4. Run GSEA
# ------------------------------

gsea_res <- GSEA(geneList     = geneList,
                 TERM2GENE    = term2gene,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

write.csv(gsea_res@result,
          "Hallmark_GSEA_results-feb13.csv",
          row.names = FALSE)

# ------------------------------
# 5. Plot
# ------------------------------

dot_p <- dotplot(gsea_res, showCategory = 20) +
  labs(title = "Hallmark GSEA: Tumor vs Adjacent") +
  theme_bw()
dot_p
ridge_plot<-ridgeplot(gsea_res, showCategory = 20)+
  labs(title = "Hallmark GSEA: Tumor vs Adjacent") +
  theme_bw()
ridge_plot

ggsave("~/work/r_project_ishu/plots/Hallmark_GSEA_dotplot-feb13.png",
       dot_p, width = 8, height = 6, dpi = 300)
table(gsea_res@result$.sign)
colnames(gsea_res@result)
table(gsea_res@result$NES > 0)
summary(gsea_res@result$NES)
gsea_res@result %>%
  filter(NES < 0) %>%
  arrange(NES)
