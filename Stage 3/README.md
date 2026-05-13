# Identify Potential Biomarkers for Cancer Detection by Analyzing Differential Expression and Machine Learning Models

**HackBio Cancer Bioinformatics Internship — Stage 3**

---

## Introduction

Lung adenocarcinoma (LUAD) is the most common histological subtype of non-small cell lung cancer (NSCLC), representing nearly half of all lung cancer cases. It originates from glandular epithelial cells lining the alveoli and is characterized by glandular structures or mucin production. Early and accurate detection remains a major clinical challenge.

This project integrates **differential expression analysis (DEA)** in R with **supervised machine learning (ML)** in Python to identify gene-level biomarkers that distinguish tumor from normal lung tissue using publicly available TCGA-LUAD data.

---

## Dataset

The Cancer Genome Atlas (TCGA) LUAD dataset was used. By June 2015, 521 LUAD samples had been analyzed and uploaded to the TCGA data portal. For this analysis, a balanced subset of **40 samples** (20 Primary Tumor + 20 Solid Tissue Normal) was selected.

| Resource | Link |
|----------|------|
| GitHub Repository | https://github.com/sanzidaanee/Hackbio-cancer-internship/tree/main/Stage%203 |
| DEG Analysis (R) | https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/Stage%203/Code/DEG.Rmd |
| ML Notebook (Python) | https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/Stage%203/Code/ML_Lung_Cancer.ipynb |
| Dataset | https://github.com/sanzidaanee/Hackbio-cancer-internship/tree/main/Stage%203/Data |
| ML Results | https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/Stage%203/Data/ml-results%20.csv |

---

## Full Project Workflow

<img width="2356" height="3256" alt="workflow_diagram (1)" src="https://github.com/user-attachments/assets/83e126ad-ac26-4013-ab04-b7c02001606a" />

The pipeline begins with the TCGA-LUAD dataset and a shared preprocessing step. It then splits into two parallel branches — the **DEG branch (R)** performs differential expression, visualization, and enrichment analysis, while the **ML branch (Python)** performs classification, feature selection, and model evaluation. Both branches converge into a combined key results section identifying the final biomarker candidates.

---

## Module 1 — Differential Expression Analysis

**Language:** R &nbsp;|&nbsp; **Tools:** `TCGAbiolinks`, `edgeR`, `biomaRt`, `gplots`  
**Full details:** [`README_DEG_Analysis.md`](./README_DEG_Analysis.md)

Identifies genes significantly differentially expressed between Primary Tumor and Solid Tissue Normal samples using a negative binomial statistical model (edgeR).

| Parameter | Value |
|-----------|-------|
| FDR-adjusted p-value | < 0.01 |
| Log₂ Fold Change (up) | > 1 |
| Log₂ Fold Change (down) | < −1 |
| Upregulated genes | **3,277** |
| Downregulated genes | **6,357** |

**Outputs:** Volcano plot · Heatmap · GO + KEGG enrichment barplots · Annotated DEG CSV files

---

## Module 2 — Machine Learning Classification

**Language:** Python &nbsp;|&nbsp; **Tools:** `scikit-learn`, `pandas`, `matplotlib`  
**Full details:** [`README_ML_Lung_Cancer.md`](./README_ML_Lung_Cancer.md)

Trains a Random Forest classifier on normalized gene expression data to classify samples and rank the top 100 discriminative genes as biomarker candidates using Recursive Feature Elimination (RFE).

| Parameter | Value |
|-----------|-------|
| Feature selection | RFE (top 100 genes) |
| Classifier | Random Forest (100 trees) |
| Train / Test split | 80% / 20% |
| Overall accuracy | **1.00** |

**Outputs:** Feature importance plot · `selected_genes.csv` · `ml-results.csv`

---

## Key Results

### Differential Expression

| Category | Count |
|----------|-------|
| Upregulated genes | 3,277 |
| Downregulated genes | 6,357 |
| Total significant DEGs | 9,634 |

### Functional Enrichment Highlights

**Upregulated in tumor:** Cell adhesion, biological adhesion, and behavior-related processes — consistent with invasion and metastatic potential in LUAD.

**Downregulated in tumor:** Cell-cell signaling, nuclear division, and mitosis — indicating disrupted cell cycle regulation and loss of normal proliferative control.

### Machine Learning

| Metric | Score |
|--------|-------|
| Overall Accuracy | 1.00 |
| Primary Tumor F1 | 1.00 |
| Solid Tissue Normal F1 | 1.00 |
| Biomarker candidates (RFE) | 100 genes |

> ⚠️ **Caveat:** Perfect accuracy on n = 40 samples may reflect overfitting. Validation on the full TCGA-LUAD cohort (521 samples) with k-fold cross-validation is strongly recommended.

---

## Conclusions

- DEG analysis successfully identified thousands of dysregulated genes in LUAD, providing insight into the molecular mechanisms driving cancer progression.
- Functional enrichment revealed biologically meaningful pathways — upregulation of adhesion/invasion and downregulation of cell cycle control — consistent with established LUAD biology.
- The Random Forest classifier with RFE identified 100 candidate biomarker genes that perfectly separated tumor from normal samples in the test set.
- Together, the two approaches provide complementary perspectives: DEA offers mechanistic insight while ML provides a data-driven ranking of the most diagnostically relevant genes.

---

## References

1. Liang Y, et al. Clinical significance of TMEM229A Q200del mutation in lung adenocarcinoma. 2023.
2. Shiba-Ishii A. Significance of stratifin in early progression of lung adenocarcinoma. *Pathology International*, 2021.
3. Herbst RS, Morgensztern D, Boshoff C. The biology and management of non-small cell lung cancer. *Nature*, 2018.
4. The Cancer Genome Atlas. Data Portal. NCI, 2015. https://tcga-data.nci.nih.gov/tcga/
5. Wang Z, Gerstein M, Snyder M. RNA-Seq: a revolutionary tool for transcriptomics. *Nature Reviews Genetics*, 2009.
6. Abdelwahab O, et al. A feature selection-based framework to identify biomarkers for cancer diagnosis. *PLOS ONE*, 2022. https://doi.org/10.1371/journal.pone.0269126

---

## Contributors

| Name | Contributions |
|------|-------------|
| Sanzida Akhter Anee | DEG and ML Analysis and report writing |

| Stéphie Raveloson | ML workflow |

| Idahosa Clinton | ML analysis |
