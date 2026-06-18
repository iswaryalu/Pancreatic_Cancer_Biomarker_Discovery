# ==========================================
# FIGURE 4 – Functional Analysis of
# DEG ∩ WGCNA Overlap Genes
# ==========================================

rm(list = ls())

# ------------------------------
# 1. Load Libraries
# ------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(dplyr)
library(pathview)
library(HDO.db)

# ------------------------------
# 2. Load Overlap Genes
# ------------------------------

overlap_file <- "~/work/r_project_ishu/Hub_DEG_overlap_genes-feb11.txt"

overlap_symbols <- read.table(overlap_file,
                              header = FALSE,
                              stringsAsFactors = FALSE)

gene_symbols <- toupper(overlap_symbols$V1)

# ------------------------------
# 3. Convert to Entrez IDs
# ------------------------------

gene_df <- bitr(gene_symbols,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

entrez_ids <- unique(gene_df$ENTREZID)

# ------------------------------
# 4. Background Universe
# (Optional but recommended)
# ------------------------------

background_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

background_df <- bitr(background_symbols,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Hs.eg.db)

background_ids <- unique(background_df$ENTREZID)

# ------------------------------
# 5. GO Enrichment (BP, CC, MF)
# ------------------------------

ego_bp <- enrichGO(gene = entrez_ids,
                   universe = background_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_cc <- enrichGO(gene = entrez_ids,
                   universe = background_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_mf <- enrichGO(gene = entrez_ids,
                   universe = background_ids,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# Save results
write.csv(as.data.frame(ego_bp),
          "~/work/r_project_ishu/KEGG_GO_COMMON_GENES/results/GO_BP_overlap.csv")

write.csv(as.data.frame(ego_cc),
          "~/work/r_project_ishu/KEGG_GO_COMMON_GENES/results/GO_CC_overlap.csv")

write.csv(as.data.frame(ego_mf),
          "~/work/r_project_ishu/KEGG_GO_COMMON_GENES/results/GO_MF_overlap.csv")

# ------------------------------
# 6. Disease Ontology (DO)
# ------------------------------

edo <- enrichDO(gene = entrez_ids,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

write.csv(as.data.frame(edo),
          "~/work/r_project_ishu/KEGG_GO_COMMON_GENES/results/DO_overlap.csv")

# ------------------------------
# 7. KEGG Enrichment
# ------------------------------

ekegg <- enrichKEGG(gene = entrez_ids,
                    organism = "hsa",
                    universe = background_ids,
                    pvalueCutoff = 0.05)

write.csv(as.data.frame(ekegg),
          "~/work/r_project_ishu/KEGG_GO_COMMON_GENES/results/KEGG_overlap.csv")

# ------------------------------
# 8. Visualization
# ------------------------------

dir.create("~/work/r_project_ishu/KEGG_GO_COMMON_GENES", showWarnings = FALSE)

# ------------------------------
# 8A. Combine BP, CC, MF
# ------------------------------

bp_df <- as.data.frame(ego_bp)
cc_df <- as.data.frame(ego_cc)
mf_df <- as.data.frame(ego_mf)

bp_df$ONTOLOGY <- "BP"
cc_df$ONTOLOGY <- "CC"
mf_df$ONTOLOGY <- "MF"

go_all <- rbind(bp_df, cc_df, mf_df)

# Take top 7 terms per ontology (like paper)
go_all <- go_all %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice(1:7) %>%
  ungroup()

# Order descriptions
go_all$Description <- factor(go_all$Description,
                             levels = rev(unique(go_all$Description)))
# ------------------------------
# 8B. Paper-style GO barplot
# ------------------------------

go_paper_plot <- ggplot(go_all,
                        aes(x = Count,
                            y = Description,
                            fill = p.adjust)) +
  geom_bar(stat = "identity") +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", switch = "y") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "GO Enrichment of Overlap Genes",
       x = "Count",
       y = NULL,
       fill = "P value") +
  theme_bw() +
  theme(
    strip.placement = "outside",   # <-- puts BP/CC/MF outside
    strip.background = element_rect(fill = "grey85"),
    strip.text.y.right = element_text(angle = 0, face = "bold"),
    panel.spacing = unit(0.5, "lines")
  )

ggsave("~/work/r_project_ishu/KEGG_GO_COMMON_GENES/plots/GO_ALL_paper_style.png",
       go_paper_plot, width = 7, height = 8)

# DO Dotplot
do_plot <- dotplot(edo, showCategory = 15) +
  ggtitle("Disease Ontology – Overlap Genes") +
  theme_bw()

ggsave("~/work/r_project_ishu/KEGG_GO_COMMON_GENES/plots/DO_dotplot.png",
       do_plot, width = 8, height = 6)

# KEGG Dotplot
kegg_plot <- dotplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment – Overlap Genes") +
  theme_bw()

ggsave("~/work/r_project_ishu/KEGG_GO_COMMON_GENES/plots/KEGG_dotplot.png",
       kegg_plot, width = 8, height = 6)
kegg_bar <- barplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment – Overlap Genes") +
  theme_bw()

ggsave("~~/work/r_project_ishu/KEGG_GO_COMMON_GENES/plots/KEGG_barplot.png",
       kegg_bar, width = 8, height = 6)

# ------------------------------
# 9. TNF Signaling Pathway Map
# ------------------------------

tnf_pathway_id <- "hsa04668"  # TNF signaling pathway

pathview(gene.data = rep(1, length(entrez_ids)),
         pathway.id = tnf_pathway_id,
         species = "hsa",
         out.suffix = "TNF_overlap")
