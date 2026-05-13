# 🫁 Machine Learning for Lung Cancer Biomarker Identification

**Project:** Identify Potential Biomarkers for Cancer Detection using ML Models  
**Dataset:** TCGA Lung Adenocarcinoma (LUAD)  
**Language:** Python  
**Notebook:** `ML_Lung_Cancer.ipynb`

---

## 📌 Overview

This module applies supervised machine learning to classify lung tissue samples as **Primary Tumor** or **Solid Tissue Normal** using gene expression data from The Cancer Genome Atlas (TCGA) LUAD dataset. By combining **Recursive Feature Elimination (RFE)** for feature selection with a **Random Forest classifier**, the pipeline identifies the most discriminating gene expression features — potential biomarkers for lung adenocarcinoma detection.

---

## 🎯 Objectives

- Predict sample type (Primary Tumor vs. Solid Tissue Normal) from gene expression profiles
- Reduce feature dimensionality using RFE to identify the most informative genes
- Train and evaluate a Random Forest classifier on the selected features
- Identify candidate biomarker genes relevant to lung cancer diagnosis

---

## 📂 Dataset

| Property | Detail |
|----------|--------|
| Source | The Cancer Genome Atlas (TCGA) |
| Cancer Type | Lung Adenocarcinoma (LUAD) |
| Total Samples | 521 (subset of 40 used: 20 Primary Tumor + 20 Solid Tissue Normal) |
| Features | Normalized gene expression values (RNA-Seq counts) |
| Target | Sample type: `Primary Tumor` / `Solid Tissue Normal` |

---

## 🔬 Workflow

```
<img width="131" height="150" alt="ml_lung_cancer_workflow" src="https://github.com/user-attachments/assets/33a08cef-a097-41e0-a291-9fab76a99439" />








Raw Gene Expression Data
        │
        ▼
Metadata Alignment (Barcode Indexing)
        │
        ▼
Common Sample Subsetting
        │
        ▼
Target Label Definition
        │
        ▼
Feature Selection (Recursive Feature Elimination - RFE)
        │
        ▼
Train/Test Split
        │
        ▼
Random Forest Classifier Training
        │
        ▼
Model Evaluation (Accuracy, Confusion Matrix, F1-Score)
        │
        ▼
Biomarker Gene Identification
```

---

## ⚙️ Methods

### 1. Data Preparation

- **Barcode indexing:** Metadata barcodes are set as the DataFrame index to ensure alignment with the gene expression matrix.
- **Sample subsetting:** Only samples present in both the metadata and gene expression dataset are retained.
- **Target definition:** The classification target is binary — Primary Tumor (1) or Solid Tissue Normal (0).

### 2. Feature Selection — Recursive Feature Elimination (RFE)

RFE iteratively removes the least important features by training a model and pruning weak predictors. This step:
- Reduces the high-dimensional gene expression space to a compact set of informative genes
- Removes noise and irrelevant features that could cause overfitting
- Produces a ranked list of genes most predictive of tumor vs. normal status

### 3. Model Training — Random Forest Classifier

A Random Forest ensemble of decision trees was trained on the RFE-selected features:

| Parameter | Value |
|-----------|-------|
| Algorithm | Random Forest |
| Task | Binary Classification |
| Input | RFE-selected gene expression features |
| Output | Primary Tumor / Solid Tissue Normal |

Random Forest is well-suited for this task because it:
- Handles high-dimensional biological data effectively
- Provides feature importance scores for biomarker ranking
- Is robust to overfitting compared to single decision trees

### 4. Model Evaluation

The classifier was evaluated using:
- **Accuracy Score**
- **Confusion Matrix**
- **Precision, Recall, and F1-Score** (per class)

---

## 📊 Results

| Metric | Primary Tumor | Solid Tissue Normal |
|--------|--------------|---------------------|
| Precision | 1.00 | 1.00 |
| Recall | 1.00 | 1.00 |
| F1-Score | 1.00 | 1.00 |
| **Overall Accuracy** | **1.00** | |

The Random Forest classifier achieved **perfect classification accuracy (100%)** on the LUAD dataset with zero misclassifications.

### Confusion Matrix

```
                     Predicted
                   Tumor | Normal
Actual  Tumor  [   TP   |   0   ]
        Normal [    0   |   TN  ]
```

> ⚠️ **Important Note on Overfitting:** The perfect accuracy should be interpreted cautiously given the small dataset size (n=40). Cross-validation on a larger cohort is recommended to validate generalizability.

---

## 🧬 Biological Significance

Genes selected by RFE as top features represent candidate biomarkers that:
- Distinguish tumor from normal lung tissue at the transcriptomic level
- May be involved in oncogenic pathways or tumor suppression
- Are consistent with findings from the parallel DEG analysis (upregulated/downregulated gene sets)

These features, when cross-referenced with differentially expressed genes, can help narrow the list of clinically relevant biomarker candidates.

---

## 📁 Output Files

| File | Description |
|------|-------------|
| `ml-results.csv` | Classification results per sample |
| `ML_Lung_Cancer.ipynb` | Full analysis notebook |

---

## 🛠️ Requirements

```bash
pip install pandas numpy scikit-learn matplotlib seaborn jupyter
```

| Package | Purpose |
|---------|---------|
| `pandas` | Data manipulation and indexing |
| `numpy` | Numerical operations |
| `scikit-learn` | RFE, Random Forest, evaluation metrics |
| `matplotlib` / `seaborn` | Visualization |

---

## ▶️ How to Run

```bash
# Clone the repository
git clone https://github.com/sanzidaanee/Hackbio-cancer-internship.git

# Navigate to the Stage 3 code directory
cd "Hackbio-cancer-internship/Stage 3/Code"

# Launch Jupyter Notebook
jupyter notebook ML_Lung_Cancer.ipynb
```

---

## 📚 References

- Abdelwahab O, Awad N, Elserafy M, Badr E (2022). A feature selection-based framework to identify biomarkers for cancer diagnosis: A focus on lung adenocarcinoma. *PLOS ONE* 17(9): e0269126.
- The Cancer Genome Atlas (2015). TCGA Data Portal. National Cancer Institute.

---

## 👥 Contributors

Sanzida Akhter Anee, Sk Arif, Nada Ghozlan, Mennatallah Mohamed Ebrahim Mahmoud, Stéphie Raveloson, Chidimma Nwaku, Zaka Ullah, Idahosa Clinton

**HackBio Cancer Bioinformatics Internship — Stage 3**
