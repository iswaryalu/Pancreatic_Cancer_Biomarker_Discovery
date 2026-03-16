# =============================
# 1. Setup
# =============================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("affy", "GEOquery", "hugene10sttranscriptcluster.db", "AnnotationDbi"), ask = FALSE)

library(affy)
library(GEOquery)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)
library(dplyr)
library(tidyverse)
library(stringr)

# =============================
# 2. Load raw CEL files
# =============================
cel_path <- "GSE62452_RAW"  # uncompressed CEL files go here
raw.data <- ReadAffy(celfile.path = cel_path)
raw.data

# =============================
# 3. RMA normalization
# =============================
normalized.data <- rma(raw.data)
normalized.expr <- as.data.frame(exprs(normalized.data))

# =============================
# 4. Map probe IDs to gene symbols
# =============================
gene_symbols <- mapIds(
  hugene10sttranscriptcluster.db,
  keys = rownames(normalized.expr),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

normalized.expr$GeneSymbol <- gene_symbols

# Remove probes without gene symbol
expr_mapped <- normalized.expr[!is.na(normalized.expr$GeneSymbol), ]

# Collapse multiple probes per gene (average)
expr_collapsed <- expr_mapped %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean))

expr_collapsed <- as.data.frame(expr_collapsed)
rownames(expr_collapsed) <- expr_collapsed$GeneSymbol
expr_collapsed <- expr_collapsed[, -1]

# Remove anything after "_" in column names
colnames(expr_collapsed) <- sapply(strsplit(colnames(expr_collapsed), "_"), `[`, 1)

# Save normalized expression matrix
write.csv(expr_collapsed, "GSE62452_RAW_normalized_expression_matrix_pancreas-tumor-val.csv")

# =============================
# 5. Extract and clean metadata from GEO
# =============================
gse <- getGEO("GSE62452", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[1]]))
head(meta)
colnames(meta)

# Clean metadata
# Clean metadata for GSE62452
meta_clean <- meta %>%
  dplyr::mutate(SampleID = geo_accession) %>%
  dplyr::select(
    SampleID,
    tissue = `tissue:ch1`,
    grade = `grading:ch1`,
    stage = `Stage:ch1`,
    survival_time_months = `survival months:ch1`,
    survival_status = `survival status:ch1`
  )
# Keep only SampleID and tissue
meta_simple <- meta_clean %>%
  dplyr::select(SampleID, tissue) %>%
  # Clean tissue values
  dplyr::mutate(tissue = case_when(
    tissue == "Pancreatic tumor" ~ "tumor",
    tissue == "adjacent pancreatic non-tumor" ~ "normal",
    TRUE ~ tissue
  ))

# Inspect
head(meta_simple)

# Save cleaned metadata
write.csv(meta_simple, "GSE62452_RAW_metadata_pancreas_tissue_only--tumor-val.csv", row.names = FALSE)
