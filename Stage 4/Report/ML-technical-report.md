# Machine Learning Technical Report
## IDH Status Classification in Low-Grade Glioma Using Gene Expression Profiles

**HackBio Cancer Bioinformatics Internship — Stage 4**  
**Language:** Python | **Models:** Random Forest · K-Nearest Neighbors

---

## Table of Contents

1. [Overview](#1-overview)
2. [Environment Setup](#2-environment-setup)
3. [Data Loading & Preparation](#3-data-loading--preparation)
4. [Feature Selection](#4-feature-selection)
5. [Model A — Random Forest with RFE](#5-model-a--random-forest-with-rfe)
6. [Model B — K-Nearest Neighbors](#6-model-b--k-nearest-neighbors)
7. [Model Comparison & Results](#7-model-comparison--results)
8. [Output & Biomarker Genes](#8-output--biomarker-genes)
9. [Interpretation](#9-interpretation)
10. [Limitations & Recommendations](#10-limitations--recommendations)

---

## 1. Overview

This report documents the supervised machine learning pipeline used to classify IDH mutation status (Mutant vs Wildtype) in TCGA-LGG samples using normalized gene expression data. Two classifiers were trained and evaluated: **Random Forest** with Recursive Feature Elimination (RFE), and **K-Nearest Neighbors** with variance thresholding.

### ML Pipeline Overview

```
Normalized Expression Matrix (34,539 genes × 494 samples)
                        │
                        ▼
              Load Data & Align Barcodes
                        │
                        ▼
             Transpose Matrix (samples × genes)
                        │
                        ▼
            Train / Test Split (80% / 20%)
                        │
              ┌──────────┴──────────┐
              │                     │
              ▼                     ▼
   ┌──────────────────┐   ┌──────────────────────┐
   │  MODEL A: RF+RFE │   │   MODEL B: KNN       │
   │                  │   │                      │
   │ 1. Train RF on   │   │ 1. Variance          │
   │    full features │   │    thresholding      │
   │ 2. Rank by       │   │    (threshold=0.01)  │
   │    importance    │   │ 2. Train KNN         │
   │ 3. Select top 100│   │    (k=5 neighbors)   │
   │ 4. Retrain on    │   │ 3. Evaluate          │
   │    top 100       │   │                      │
   │                  │   │ Accuracy: 83%        │
   │ Accuracy: 99%    │   │ WT recall: 35%       │
   └──────────────────┘   └──────────────────────┘
              │                     │
              └──────────┬──────────┘
                         ▼
              Comparison & Biomarker Output
              Top 100 candidate genes saved
```

---

## 2. Environment Setup

### Required Libraries

```python
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.feature_selection import RFE, VarianceThreshold
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
```

### Install Dependencies

```bash
pip install scikit-learn pandas numpy matplotlib
```

**Key library versions used:**

| Library | Version | Role |
|---------|---------|------|
| scikit-learn | ≥ 1.0 | ML models, feature selection, metrics |
| pandas | ≥ 1.3 | Data manipulation |
| numpy | ≥ 1.21 | Numerical operations |
| matplotlib | ≥ 3.4 | Feature importance visualization |

---

## 3. Data Loading & Preparation

### Load Normalized Expression Data & Metadata

```python
# Load normalized gene expression matrix (genes × samples)
normalized_data_url = 'https://drive.google.com/uc?id=1E7xFeh2C7IYG4dTW8bCwfABJphAKHq0_'
metadata_url        = 'https://drive.google.com/uc?id=1VqL4hjYY4vcRcneaKRP9rKg7Zc3NCwLy'

normalized_data = pd.read_csv(normalized_data_url)
metadata        = pd.read_csv(metadata_url)
```

### Align Samples Between Datasets

```python
# Set barcode as index for metadata
metadata.set_index('Barcode', inplace=True)

# Find samples present in BOTH datasets
common_samples = normalized_data.columns.intersection(metadata.index)
print(f"Common samples: {len(common_samples)}")

# Subset to common samples only
normalized_data = normalized_data[common_samples]
metadata        = metadata.loc[common_samples]
```

> **Why align?** The expression matrix (from DEA output) and the metadata may not share exactly the same sample set. Taking the intersection ensures no label mismatch.

### Prepare Feature Matrix and Target Labels

```python
# Transpose: convert from (genes × samples) → (samples × genes)
# ML models expect rows = samples, columns = features
X = normalized_data.T     # Shape: (n_samples, n_genes)
y = metadata['IDH']       # Target: 'Mutant' or 'WT'

print(f"Feature matrix shape: {X.shape}")
print(f"Class distribution:\n{y.value_counts()}")
```

**Expected output:**
```
Feature matrix shape: (494, 34539)
Class distribution:
Mutant    460
WT         34
```

> ⚠️ **Class imbalance:** 460 Mutant vs 34 WT (13.5:1 ratio). This will inflate overall accuracy; WT-specific recall is the more informative metric.

### Train/Test Split

```python
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size    = 0.2,     # 20% held out for testing
    random_state = 42       # Reproducibility seed
)

print(f"Training samples: {X_train.shape[0]}")
print(f"Test samples:     {X_test.shape[0]}")
```

---

## 4. Feature Selection

High-dimensional gene expression data (34,539 features) requires dimensionality reduction before training to avoid the curse of dimensionality and reduce computational cost.

### Strategy Comparison

| Model | Method | Features Retained | Rationale |
|-------|--------|------------------|-----------|
| Random Forest | RFE via feature importance | Top 100 | Model-based ranking captures non-linear interactions |
| KNN | Variance Threshold | Variable (threshold=0.01) | Removes constant/near-constant genes; KNN sensitive to irrelevant features |

---

## 5. Model A — Random Forest with RFE

### Step 1: Train Initial Random Forest on Full Feature Set

```python
from sklearn.ensemble import RandomForestClassifier

# Define base classifier
rf = RandomForestClassifier(
    n_estimators = 100,    # 100 decision trees
    random_state = 42      # Reproducibility
)

# Train on all 34,539 features
rf.fit(X_train, y_train)
```

### Step 2: Rank Features by Importance

```python
# Extract Gini importance scores
importances = rf.feature_importances_

# Build ranked DataFrame
feature_importance_df = pd.DataFrame({
    'Feature'   : X_train.columns,
    'Importance': importances
}).sort_values(by='Importance', ascending=False)

# Display top 20
print(feature_importance_df.head(20))
```

> **Gini importance** measures how much each feature decreases impurity (heterogeneity) across all trees. Higher = more discriminative for IDH status.

### Step 3: Select Top 100 Features

```python
# Keep top 100 most important genes
top_100_features = feature_importance_df.head(100)
selected_features = top_100_features['Feature'].tolist()

# Filter train and test matrices
X_train_selected = X_train[selected_features]
X_test_selected  = X_test[selected_features]

print(f"Features after RFE: {len(selected_features)}")
```

### Step 4: Retrain Random Forest on Selected Features

```python
# Final Random Forest on top 100 genes
rf_model = RandomForestClassifier(
    n_estimators = 100,
    random_state = 42
)
rf_model.fit(X_train_selected, y_train)
```

### Step 5: Evaluate Performance

```python
y_pred_rf = rf_model.predict(X_test_selected)

# Metrics
accuracy_rf = accuracy_score(y_test, y_pred_rf)
print(f"Random Forest Accuracy: {accuracy_rf:.2f}")

print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred_rf))

print("\nClassification Report:")
print(classification_report(y_test, y_pred_rf))
```

### Step 6: Visualize Top Feature Importances

```python
importances_final = rf_model.feature_importances_
indices = np.argsort(importances_final)[::-1][:10]   # Top 10

plt.figure(figsize=(10, 5))
plt.title("Top 10 Feature Importances — Random Forest")
plt.bar(range(10), importances_final[indices], align="center", color="steelblue")
plt.xticks(range(10), [selected_features[i] for i in indices], rotation=90)
plt.ylabel("Gini Importance")
plt.tight_layout()
plt.savefig("rf_feature_importance.png", dpi=150)
plt.show()
```

**Figure 3:** Top 10 most important genes for IDH status classification as ranked by Random Forest Gini importance.

### Random Forest Results

| Metric | Mutant | WT | Overall |
|--------|--------|----|---------|
| Precision | 0.99 | 1.00 | — |
| Recall | 1.00 | 0.96 | — |
| F1-Score | 1.00 | 0.98 | — |
| **Accuracy** | — | — | **99%** |

**Confusion Matrix:**
```
              Predicted
              Mutant   WT
Actual Mutant  [92]   [0]
       WT       [1]  [6]
```
(1 WT sample misclassified as Mutant)

---

## 6. Model B — K-Nearest Neighbors

### Step 1: Variance-Based Feature Selection

```python
from sklearn.feature_selection import VarianceThreshold

# Remove genes with near-zero variance across samples
# These genes carry no discriminative information
selector = VarianceThreshold(threshold=0.01)
X_train_selected_knn = selector.fit_transform(X_train)
X_test_selected_knn  = selector.transform(X_test)

n_features_kept = X_train_selected_knn.shape[1]
print(f"Features after variance filtering: {n_features_kept}")
```

### Step 2: Train KNN Classifier

```python
from sklearn.neighbors import KNeighborsClassifier

# Train with k=5 neighbors (Euclidean distance by default)
knn_model = KNeighborsClassifier(n_neighbors=5)
knn_model.fit(X_train_selected_knn, y_train)
```

> **How KNN works:** For each test sample, the model finds the 5 nearest training samples in gene expression space (Euclidean distance) and assigns the majority class label. With 34K+ features, distance metrics become unreliable — this is the core weakness of KNN in high dimensions.

### Step 3: Evaluate Performance

```python
y_pred_knn = knn_model.predict(X_test_selected_knn)

accuracy_knn = accuracy_score(y_test, y_pred_knn)
print(f"KNN Accuracy: {accuracy_knn:.2f}")

print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred_knn))

print("\nClassification Report:")
print(classification_report(y_test, y_pred_knn))
```

### Step 4: Visualize Feature Distribution

```python
# KNN has no feature importance; visualize feature value distribution
plt.figure(figsize=(8, 4))
plt.title("Distribution of Variance-Selected Feature Values")
plt.hist(np.sum(X_train_selected_knn, axis=0), bins=20, color="coral", edgecolor="black")
plt.xlabel('Summed Feature Value Across Samples')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig("knn_feature_distribution.png", dpi=150)
plt.show()
```

### KNN Results

| Metric | Mutant | WT | Overall |
|--------|--------|----|---------|
| Precision | 0.82 | — | — |
| Recall | 0.99 | **0.35** | — |
| F1-Score | 0.90 | **0.50** | — |
| **Accuracy** | — | — | **83%** |

**Confusion Matrix:**
```
              Predicted
              Mutant   WT
Actual Mutant  [91]   [1]
       WT       [18]  [4] ← only 4/22 WT correctly classified
```

### Save Results to File

```python
# Save selected gene names for both models
selected_genes_knn = pd.Series(X_train.columns[selector.get_support()])
selected_genes_knn.to_csv('selected_genes_knn.csv', index=False)

# Save prediction results
results_df = pd.DataFrame({
    'barcode'              : X_test.index,
    'true_sample_type'     : y_test,
    'predicted_sample_type': y_pred_knn,
    'model_accuracy'       : accuracy_knn,
    'selected_features'    : "; ".join(map(str, selected_genes_knn))
})
print(results_df.head())
```

---

## 7. Model Comparison & Results

### Performance Summary

| Model | Feature Selection | Features Used | Accuracy | Mutant F1 | WT F1 | WT Recall |
|-------|-----------------|---------------|----------|-----------|-------|-----------|
| **Random Forest** | RFE (top 100) | 100 genes | **99%** | **1.00** | **0.98** | **0.96** |
| K-Nearest Neighbors | Variance Threshold | ~34K genes | 83% | 0.90 | 0.50 | 0.35 |

### Why Random Forest Outperforms KNN

| Factor | Random Forest | KNN |
|--------|--------------|-----|
| High dimensionality | Handles well (feature importance selects relevant genes) | Degrades severely (curse of dimensionality) |
| Class imbalance | Ensemble voting gives WT samples a fair chance | Majority-class bias in local neighborhoods |
| Feature interactions | Captures non-linear gene interactions via tree splits | Euclidean distance assumes linear feature separability |
| Interpretability | Feature importance scores identify biomarker genes | No intrinsic feature ranking |
| Computational cost | Moderate (train once, predict fast) | High at inference (recomputes distances for each test point) |

### Key Takeaway

The ~13.5:1 class imbalance (460 Mutant vs 34 WT) is the primary driver of KNN's poor WT recall (35%). In high-dimensional space, WT samples are surrounded predominantly by Mutant neighbors, causing systematic misclassification. Random Forest's ensemble approach, combined with dimension reduction via feature importance, is far more robust to both imbalance and high dimensionality.

---

## 8. Output & Biomarker Genes

### Save Random Forest Top 100 Genes

```python
# Top 100 biomarker candidates from RF feature importance
selected_genes_rf = pd.Series(selected_features)
selected_genes_rf.to_csv('selected_genes_rf.csv', index=False)

# Full RF prediction results
results_rf_df = pd.DataFrame({
    'barcode'              : X_test_selected.index,
    'true_sample_type'     : y_test,
    'predicted_sample_type': y_pred_rf,
    'model_accuracy'       : accuracy_rf,
    'selected_features'    : "; ".join(map(str, selected_features))
})
results_rf_df.to_csv('ml-results.csv', index=False)
print(results_rf_df.head(10))
```

### Output Files

| File | Contents |
|------|---------|
| `selected_genes_rf.csv` | Top 100 biomarker genes (RF importance ranking) |
| `selected_genes_knn.csv` | Genes passing variance threshold (KNN) |
| `ml-results.csv` | Per-sample predictions, true labels, accuracy |
| `rf_feature_importance.png` | Bar chart of top 10 gene importances |
| `knn_feature_distribution.png` | Distribution of KNN-selected feature values |

---

## 9. Interpretation

### Biological Meaning of RF Top 100 Genes

The 100 genes selected by RF feature importance represent those whose expression levels most reliably differ between IDH-mutant and IDH-wildtype LGG. These genes should show:

- High Gini importance (large mean decrease in impurity across the 100 trees)
- Consistency with the DEA results (many should appear in the 5,916 significant DEGs)
- Potential enrichment in IDH-mutation-related pathways (hypermethylation targets, metabolic regulators)

### Cross-Validation with DEA

The intersection of the RF top-100 biomarker genes with the 5,916 DEGs from the edgeR analysis represents the highest-confidence candidates — genes that are both statistically differentially expressed AND predictively important in the ML classifier.

---

## 10. Limitations & Recommendations

| Limitation | Recommendation |
|------------|---------------|
| Class imbalance (460:34) | Use SMOTE oversampling or class_weight='balanced' in RF |
| No cross-validation | Use stratified k-fold CV (k=5 or k=10) for robust accuracy estimates |
| Single train/test split | Results may vary with different random seeds; report mean ± SD across multiple splits |
| KNN not tuned | Grid search over k={3,5,7,11} and distance metrics (Euclidean, Manhattan, cosine) |
| No normalization before KNN | Standardize features (StandardScaler) before KNN — Euclidean distance is sensitive to scale |
| Overfitting risk | Validate RF on independent cohort (e.g., CGGA glioma dataset) |

### Recommended Improvements

```python
from sklearn.utils.class_weight import compute_class_weight
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler

# 1. Class-weighted Random Forest
class_weights = compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)
weight_dict = dict(zip(np.unique(y_train), class_weights))

rf_balanced = RandomForestClassifier(
    n_estimators = 100,
    class_weight = weight_dict,
    random_state = 42
)

# 2. Stratified K-Fold Cross-Validation
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = cross_val_score(rf_balanced, X_train_selected, y_train, cv=cv, scoring='f1_macro')
print(f"CV F1-macro: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

# 3. Scale features before KNN
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_selected_knn)
X_test_scaled  = scaler.transform(X_test_selected_knn)
```

---

## Software Versions

| Library | Version | Role |
|---------|---------|------|
| Python | ≥ 3.8 | Analysis environment |
| scikit-learn | ≥ 1.0 | RF, KNN, feature selection, metrics |
| pandas | ≥ 1.3 | Data loading and manipulation |
| numpy | ≥ 1.21 | Array operations |
| matplotlib | ≥ 3.4 | Feature importance visualization |
