# =============================
# 1. Setup
# =============================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install necessary Bioconductor packages
BiocManager::install(c("affy", "GEOquery", "hgu133plus2.db", "AnnotationDbi"), ask = FALSE)

# Load libraries
library(affy)
library(GEOquery)
library(hgu133plus2.db)  # HG-U133 Plus 2.0 annotation
library(AnnotationDbi)
library(dplyr)
library(tidyverse)

# =============================
# 2. Load raw CEL files
# =============================

cel_path <- "GSE15932_RAW"  # Set this to your folder with CEL files
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
  hgu133plus2.db,
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
# Remove the ".CEL.gz" suffix from column names
colnames(expr_collapsed) <- sub("\\.CEL\\.gz$", "", colnames(expr_collapsed))

# Check the result
colnames(expr_collapsed)

# Clean sample names
colnames(expr_collapsed) <- sapply(strsplit(colnames(expr_collapsed), "_"), `[`, 1)

# Save normalized expression matrix
write.csv(expr_collapsed, "normalized_expression_matrix_GSE15932.csv")

# =============================
# 5. Extract and clean metadata from GEO
# =============================
gse <- getGEO("GSE15932", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[1]]))

# Inspect metadata columns
head(meta)
colnames(meta)
# Create a Group column based on title
meta_clean <- meta %>%
  select(geo_accession, title, `diagnosis:ch1`, `sex:ch1`, `age of patient:ch1`) %>%
  rename(SampleID = geo_accession,
         Diagnosis = `diagnosis:ch1`,
         Sex = `sex:ch1`,
         Age = `age of patient:ch1`) %>%
  mutate(Group = case_when(
    grepl("pancreatic cancer and diabetes", title, ignore.case = TRUE) ~ "PC+Diabetes",
    grepl("diabetes mellitus", title, ignore.case = TRUE) &
      !grepl("pancreatic cancer", title, ignore.case = TRUE) ~ "Diabetes only",
    grepl("pancreatic cancer", title, ignore.case = TRUE) &
      grepl("without diabetes", title, ignore.case = TRUE) ~ "PC only",
    grepl("healthy control", title, ignore.case = TRUE) ~ "Healthy",
    TRUE ~ "Unknown"
  ))

# Check the result
table(meta_clean$Group)


# Check the result
table(meta_clean$Group)

# Save cleaned metadata
write.csv(meta_clean, "metadata_GSE15932_clean.csv", row.names = FALSE)

