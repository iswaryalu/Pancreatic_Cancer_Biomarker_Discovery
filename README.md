# Pancreatic Cancer Biomarker Discovery

## Integrative Co-expression Network and Machine Learning Analysis for Early Non-Invasive Detection of Pancreatic Cancer

### Project Summary

Pancreatic cancer is one of the most lethal malignancies worldwide, largely due to late-stage diagnosis and limited availability of reliable early detection biomarkers. This project aimed to identify robust transcriptomic biomarkers associated with pancreatic cancer by integrating differential gene expression analysis, weighted gene co-expression network analysis (WGCNA), functional enrichment, and machine learning approaches.

Using publicly available microarray datasets from the Gene Expression Omnibus (GEO), biologically relevant genes were identified, prioritized through multiple machine learning algorithms, and validated across independent tissue and blood-based datasets. The study demonstrates the potential of combining network biology and machine learning for non-invasive cancer biomarker discovery.

---

## Objectives

* Identify differentially expressed genes associated with pancreatic cancer.
* Construct gene co-expression networks to identify disease-related modules and hub genes.
* Integrate differential expression and network analysis results to obtain biologically significant candidate genes.
* Perform functional enrichment and pathway analysis to understand disease mechanisms.
* Apply machine learning techniques for feature selection and biomarker prioritization.
* Develop and validate a diagnostic model for early pancreatic cancer detection.

---

## Datasets

The analysis was performed using publicly available Affymetrix microarray datasets obtained from the NCBI Gene Expression Omnibus (GEO).

| Dataset   | Purpose                                                           |
| --------- | ----------------------------------------------------------------- |
| GSE183795 | Discovery cohort for differential expression and network analysis |
| GSE62452  | Independent tumour tissue validation                              |
| GSE15932  | Peripheral blood validation cohort                                |

---

## Methodology

### 1. Data Preprocessing

* Downloaded Affymetrix microarray datasets from GEO.
* Performed quality control and Robust Multi-array Average (RMA) normalization.
* Mapped probe identifiers to gene symbols.
* Generated normalized expression matrices for downstream analysis.

### 2. Differential Gene Expression Analysis

* Conducted using the Limma package in R.
* Applied statistical filtering based on adjusted p-values and fold-change thresholds.
* Identified significantly upregulated and downregulated genes.
* Visualized results using volcano plots and heatmaps.

### 3. Weighted Gene Co-expression Network Analysis (WGCNA)

* Constructed a signed gene co-expression network.
* Identified co-expressed gene modules associated with pancreatic cancer.
* Determined module-trait relationships.
* Extracted hub genes from disease-associated modules.

### 4. Functional Enrichment Analysis

* Gene Ontology (GO) enrichment analysis.
* KEGG pathway enrichment analysis.
* Biological interpretation of candidate genes and pathways involved in pancreatic tumorigenesis.

### 5. Machine Learning-Based Biomarker Selection

Multiple feature selection algorithms were employed:

* LASSO Regression
* Support Vector Machine Recursive Feature Elimination (SVM-RFE)
* Random Forest

Genes consistently identified across algorithms were selected as high-confidence biomarkers.

### 6. Model Development and Validation

* Developed a logistic regression diagnostic model.
* Evaluated performance using Receiver Operating Characteristic (ROC) analysis.
* Validated findings using independent tissue and blood datasets.

---

## Key Results

### Differential Expression Analysis

* Identified **311 Differentially Expressed Genes (DEGs)**

  * 193 Upregulated Genes
  * 118 Downregulated Genes

### Network Analysis

* Identified **225 Hub Genes** from disease-associated co-expression modules.
* Significant tumour-associated modules were prioritized for downstream analysis.

### Integrated Candidate Genes

* Obtained **79 overlapping genes** through integration of DEGs and WGCNA hub genes.

### Final Biomarker Signature

The following five genes were consistently identified through machine learning-based feature selection:

* AHNAK2
* TRIM29
* SLC6A14
* ITGB4
* CTSE

---

## Model Performance

| Validation Dataset                  | AUC   |
| ----------------------------------- | ----- |
| Independent Tumour Cohort           | 0.93  |
| Peripheral Blood Cohort             | 0.80  |
| Blood Cohort (Non-diabetic Samples) | ~0.91 |

The model demonstrated strong diagnostic performance and supports the feasibility of translating tumour-derived transcriptomic signatures into blood-based diagnostic applications.

---

## Bioinformatics Workflow

GEO Datasets
→ Data Preprocessing & Normalization
→ Differential Expression Analysis (Limma)
→ WGCNA
→ Functional Enrichment Analysis
→ Machine Learning Feature Selection
→ Biomarker Discovery
→ Logistic Regression Model Development
→ External Validation

---

## Tools & Technologies

### Programming & Analysis

* R
* Bioconductor

### Bioinformatics Packages

* GEOquery
* Limma
* WGCNA
* clusterProfiler
* AnnotationDbi
* STRINGdb

### Machine Learning

* glmnet (LASSO)
* e1071 (SVM)
* randomForest
* pROC

### Visualization

* ggplot2
* pheatmap
* EnhancedVolcano

---

## Repository Structure

```text
DEG_and_WGCNA_analysis/
│
├── deg/
├── wgcna/
│
enrichment_analysis/
│
ml/
├── lasso/
├── svm_rfe/
├── random_forest/
│
results/
figures/
```

---

## Skills Demonstrated

* Transcriptomics Data Analysis
* Gene Expression Profiling
* Differential Expression Analysis
* Weighted Gene Co-expression Network Analysis (WGCNA)
* Functional Enrichment Analysis
* Biomarker Discovery
* Machine Learning for Biological Data
* Statistical Data Analysis
* Data Visualization
* Reproducible Bioinformatics Workflows

---

## Author

**Iswarya L U**

M.Tech Computational Biology
Pondicherry University

Bioinformatics | Computational Biology | Transcriptomics | Functional Genomics | Machine Learning for Life Sciences
