# Differential Gene Expression Profiling and Unsupervised Clustering to Predict Glioma Prognosis

**HackBio Cancer Bioinformatics Internship — Stage 4**

---

## Project Overview

Gliomas are primary brain tumors arising from glial cells, accounting for approximately 80% of all malignant brain tumors. This project integrates **differential expression analysis (DEA)** in R with **supervised machine learning (ML)** in Python to characterize IDH mutation-driven transcriptomic differences in low-grade glioma (LGG) and build predictive models that classify IDH status from gene expression data.

**Key biological question:** Can transcriptomic profiles reliably distinguish IDH-mutant from IDH-wildtype low-grade gliomas, and which genes serve as the strongest predictive biomarkers?

---

## Dataset

| Property | Details |
|----------|---------|
| Source | The Cancer Genome Atlas (TCGA) |
| Project | TCGA-LGG (Low-Grade Glioma) |
| Total samples | 534 patients |
| IDH Mutant | 460 samples |
| IDH Wildtype | 34 samples |
| Total genes | 60,660 |
| Filtered genes used | 34,539 |
| Data type | RNA-Seq (STAR Counts) |

---

## Full Project Workflow

<img width="3276" height="4716" alt="workflow_diagram" src="https://github.com/user-attachments/assets/ec61153d-1a53-4e62-b576-40f3bed56126" />

The pipeline begins with the TCGA-LGG dataset and a shared preprocessing step. It then splits into two parallel branches — the **DEA branch (R)** performs differential expression, visualization, and GO enrichment analysis, while the **ML branch (Python)** performs feature selection, classification, and model evaluation. Both branches converge into a combined key results section identifying the final biomarker candidates.

---

## Module 1 — Differential Expression Analysis (R)

**Tools:** `TCGAbiolinks` · `edgeR` · `EDASeq` · `biomaRt` · `gplots`  
**Full details:** [`DEA_Technical_Report.md`](./DEA_Technical_Report.md)

### Analysis Parameters

| Parameter | Value |
|-----------|-------|
| Normalization method | Gene-length (quantile) |
| FDR cutoff | < 0.01 |
| log₂ Fold Change (up) | > 1 |
| log₂ Fold Change (down) | < −1 |
| DE pipeline | edgeR (negative binomial) |

### Results Summary

| Category | Count |
|----------|-------|
| Filtered genes analyzed | 34,539 |
| Total significant DEGs | **5,916** |
| Upregulated (Mutant vs WT) | **1,681** |
| Downregulated (Mutant vs WT) | **4,235** |

### Key Pathway Enrichments

**Upregulated in IDH-Mutant:**
- Synaptic transmission
- Cell-cell signaling
- Homophilic cell adhesion
- Transmission of nerve impulse

**Downregulated in IDH-Mutant:**
- Anterior/posterior pattern formation
- Skeletal system development
- Embryonic morphogenesis
- Regionalization & pattern specification

---

## Module 2 — Machine Learning Classification (Python)

**Tools:** `scikit-learn` · `pandas` · `numpy` · `matplotlib`  
**Full details:** [`ML_Technical_Report.md`](./ML_Technical_Report.md)

### Models Evaluated

| Model | Feature Selection | Accuracy | Mutant F1 | WT F1 |
|-------|------------------|----------|-----------|-------|
| Random Forest | RFE (top 100 genes) | **99%** | 1.00 | 0.98 |
| K-Nearest Neighbors | Variance Threshold | **83%** | 0.90 | 0.50 |

### Key Findings

- **Random Forest** vastly outperformed KNN, especially for wildtype classification
- RFE selected 100 high-importance biomarker genes from the full expression matrix
- KNN recall for WT samples was only 35%, indicating poor sensitivity for the minority class
- RF misclassified only 1 sample across the entire test set

> ⚠️ **Caveat:** Near-perfect RF accuracy on n=534 samples with strong class imbalance (460 Mutant vs 34 WT) warrants validation with cross-validation and independent cohorts.

---

## Key Conclusions

1. **IDH mutation status** is the primary molecular stratifier in LGG, producing a highly distinct transcriptomic fingerprint
2. **Downregulated genes** in IDH-mutant tumors are enriched in developmental/morphogenetic pathways — consistent with epigenetic silencing via hypermethylation
3. **Upregulated genes** include neural signaling and adhesion pathways, potentially reflecting altered tumor-microenvironment interactions
4. **Random Forest with RFE** is the superior classifier, achieving 99% accuracy and identifying 100 candidate diagnostic biomarkers
5. **Class imbalance** (460 Mutant vs 34 WT) inflates accuracy metrics — weighted evaluation metrics and resampling strategies are recommended for future work

---

## References

1. Ahir BK, Engelhard HH, Lakka SS. Tumor development and angiogenesis in adult brain tumors: glioblastoma. *Molecular Neurobiology* 2020; 57: 2461–2478.
2. Nakhate V, Lasica AB, Wen PY. The Role of Mutant IDH Inhibitors in the Treatment of Glioma. *Current Neurology and Neuroscience Reports* 2024.
3. The Cancer Genome Atlas. Data Portal. NCI, 2015. https://tcga-data.nci.nih.gov/tcga/
4. McCarthy DJ, Chen Y, Smyth GK. Differential expression analysis of multifactor RNA-Seq experiments. *Nucleic Acids Research* 2012; 40(10): 4288–4297.
5. Colaprico A et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic Acids Research* 2016; 44(8): e71.
6. Liaw A, Wiener M. Classification and Regression by randomForest. *R News* 2002; 2(3): 18–22.

---

## Contributors

| Name | Slack Handle |
|------|-------------|
| Sanzida Akhter Anee | @Sanzida |
| Nada Ghozlan | @Nad1 |
| Stéphie Raveloson | @StephieRav |
| Chidimma Nwaku | @Mma |
| Idahosa Clinton | @doc_idahosa |
