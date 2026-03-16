# =============================
# 1. Setup
# =============================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("affy", "GEOquery", "hugene10sttranscriptcluster.db", "AnnotationDbi"), ask = FALSE)

# Load libraries
library(affy)
library(GEOquery)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)
library(dplyr)
library(tidyverse)

# =============================
# 2. Load raw CEL files
# =============================
cel_path <- "GSE49641_RAW"
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

# Add gene symbols to expression matrix
normalized.expr$GeneSymbol <- gene_symbols

# Remove probes without a gene symbol
expr_mapped <- normalized.expr[!is.na(normalized.expr$GeneSymbol), ]

# Collapse multiple probes per gene (average)
expr_collapsed <- expr_mapped %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean))

# Set GeneSymbol as rownames and remove column
expr_collapsed <- as.data.frame(expr_collapsed)
rownames(expr_collapsed) <- expr_collapsed$GeneSymbol
expr_collapsed <- expr_collapsed[, -1]
colnames(expr_collapsed) <- sapply(strsplit(colnames(expr_collapsed), "_"), `[`, 1)
# Save normalized expression matrix
write.csv(expr_collapsed, "normalized_expression_matrix_blood.csv")

# =============================
# 5. Extract and clean metadata from GEO
# =============================
gse <- getGEO("GSE49641", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[1]]))

# Inspect metadata columns to see which contain relevant info
head(meta)
colnames(meta)

# =============================
# 1. Clean metadata for GSE49641
# =============================

# Make a copy of the metadata

library(dplyr)

# Convert to tibble first
meta_clean <- as_tibble(meta)

library(dplyr)

meta_clean <- meta %>%
  dplyr::mutate(SampleID = sapply(strsplit(geo_accession, "_"), `[`, 1)) %>%
  dplyr::select(
    SampleID,
    age = `age:ch1`,
    gender = `gender:ch1`,
    cell_type = `cell type:ch1`,
    disease_state = `disease state:ch1`,
    pathological_stage = `pathological staging:ch1`
  )

# Clean column names
colnames(meta_clean) <- colnames(meta_clean) %>%
  stringr::str_replace_all(" ", "_") %>%
  stringr::str_replace_all(":", "_") %>%
  stringr::str_replace_all("-", "_")

head(meta_clean)


# Save cleaned metadata
write.csv(meta_clean, "metadata_blood_clean.csv", row.names = FALSE)

