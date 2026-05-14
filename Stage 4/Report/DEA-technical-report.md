# Differential Expression Analysis — Technical Report
## Glioma IDH Status Classification from TCGA-LGG RNA-Seq Data

**HackBio Cancer Bioinformatics Internship — Stage 4**  
**Language:** R | **Pipeline:** TCGAbiolinks + edgeR

---

## Table of Contents

1. [Overview](#1-overview)
2. [Environment Setup](#2-environment-setup)
3. [Data Acquisition](#3-data-acquisition)
4. [Metadata Construction & IDH Status](#4-metadata-construction--idh-status)
5. [Normalization & Filtering](#5-normalization--filtering)
6. [Differential Expression Analysis](#6-differential-expression-analysis)
7. [Visualization](#7-visualization)
8. [Functional Enrichment Analysis](#8-functional-enrichment-analysis)
9. [Results Summary](#9-results-summary)
10. [Interpretation](#10-interpretation)

---

## 1. Overview

This analysis characterizes transcriptomic differences between IDH-mutant and IDH-wildtype low-grade gliomas (LGG) using RNA-Seq data from the TCGA-LGG cohort. The edgeR negative binomial model is applied after gene-length normalization to identify differentially expressed genes (DEGs), followed by GO biological process enrichment analysis.

### Analysis Pipeline

```
TCGA-LGG Download
       │
       ▼
Metadata Construction
(IDH status extraction: Mutant / WT)
       │
       ▼
Raw Count Matrix (unstranded)
       │
       ▼
Gene-Length Normalization (TCGAanalyze_Normalization)
       │
       ▼
Quantile Filtering (bottom 25% removed)
       │
       ▼
Differential Expression Analysis (edgeR)
FDR < 0.01, |logFC| > 1
       │
       ▼
DEG Classification (UP / DOWN / NS)
       │
       ├──► Volcano Plot
       ├──► Heatmap (hierarchical clustering)
       └──► GO Enrichment (Biological Process)
```

---

## 2. Environment Setup

### Install Required Packages

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")   # TCGA data access & DE analysis
BiocManager::install("edgeR")          # Differential expression
BiocManager::install("EDASeq")         # Normalization
BiocManager::install("SummarizedExperiment")  # Data structure
BiocManager::install("biomaRt")        # Gene annotation

# CRAN packages
install.packages("dplyr")
install.packages("gplots")
```

### Load Libraries

```r
library(TCGAbiolinks)         # Accessing and querying TCGA data
library(edgeR)                # Differential expression analysis
library(EDASeq)               # Exploratory data analysis and normalization
library(SummarizedExperiment) # Access assay data
library(biomaRt)              # Access BioMart databases
library(dplyr)                # Data manipulation
library(gplots)               # Heatmap generation
```

---

## 3. Data Acquisition

### Query TCGA-LGG RNA-Seq Data

```r
# List all GDC projects (for reference)
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-LGG')

# Build query for LGG gene expression data
lggQ <- GDCquery(
  project             = 'TCGA-LGG',
  data.category       = 'Transcriptome Profiling',
  experimental.strategy = "RNA-Seq",
  workflow.type       = "STAR - Counts",
  access              = "open",
  data.type           = "Gene Expression Quantification"
)

# Download and prepare
GDCdownload(lggQ)
lgg.data <- GDCprepare(lggQ)
```

### Download IDH Mutation Data

```r
mutation_query <- GDCquery(
  project       = "TCGA-LGG",
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(mutation_query)
mutation_data <- GDCprepare(mutation_query)

# Filter for IDH1 and IDH2 only
idh_mutations <- mutation_data[mutation_data$Hugo_Symbol %in% c("IDH1", "IDH2"), ]
```

**Dataset dimensions after download:**
- Genes: 60,660
- Samples: 534 (460 Mutant, 34 WT, 40 NA/excluded)

---

## 4. Metadata Construction & IDH Status

### Build Metadata DataFrame

```r
Metadata <- data.frame(
  "Barcode"              = lgg.data$barcode,
  "Tumor_Sample_Barcode" = lgg.data$bcr_patient_barcode,
  "tumor_type"           = lgg.data$tumor_descriptor,
  "sample_id"            = lgg.data$sample_id,
  "sample"               = lgg.data$sample_type,
  "gender"               = lgg.data$gender,
  "IDH"                  = lgg.data$paper_IDH.status,        # KEY: Mutant / WT
  "Mutation"             = lgg.data$paper_Mutation.Count,
  "TERT"                 = lgg.data$paper_TERT.expression.status,
  "TERT_status"          = lgg.data$paper_TERT.promoter.status,
  "ATRX_status"          = lgg.data$paper_ATRX.status,
  "RNA_cluster"          = lgg.data$paper_IDH.specific.RNA.Expression.Cluster,
  "Random_Forest_cluster"= lgg.data$paper_Random.Forest.Sturm.Cluster
)

# Save metadata
write.csv(Metadata, "TCGA_LGG_metadata.csv", row.names = FALSE)
```

### Align Mutation Data with RNA-Seq Metadata

```r
# Standardize patient IDs to first 12 characters
idh_mutations$Tumor_Sample_Barcode <- substr(idh_mutations$Tumor_Sample_Barcode, 1, 12)
Metadata$Tumor_Sample_Barcode <- substr(Metadata$Tumor_Sample_Barcode, 1, 12)

# Merge datasets
merged_data <- merge(Metadata, idh_mutations, by = "Tumor_Sample_Barcode")

# Separate by IDH gene
idh1_data <- merged_data[merged_data$Hugo_Symbol == "IDH1", ]
idh2_data <- merged_data[merged_data$Hugo_Symbol == "IDH2", ]

# Separate by mutation status
mutant_samples  <- merged_data[merged_data$IDH == "Mutant", ]
wildtype_samples <- merged_data[merged_data$IDH == "WT", ]
```

---

## 5. Normalization & Filtering

### Extract Raw Count Matrix

```r
lgg.raw.data <- assays(lgg.data)
dim(lgg.raw.data$unstranded)   # 60660 × 534

# Select samples with known IDH status (WT + Mutant only)
selectedBarcodes <- c(
  subset(Metadata, IDH == "WT")$Barcode,
  subset(Metadata, IDH == "Mutant")$Barcode
)

selectedData <- lgg.raw.data$unstranded[, selectedBarcodes]
dim(selectedData)   # 60660 × (34 + 460) = 60660 × 494
```

### Gene-Length Normalization

```r
# Load gene length information bundled with TCGAbiolinks
data(geneInfoHT, package = "TCGAbiolinks")

# Normalize by gene length (RPKM-like normalization)
normData <- TCGAanalyze_Normalization(
  tabDF    = selectedData,
  geneInfo = geneInfoHT,
  method   = "geneLength"
)

# Save normalized data for downstream ML use
write.csv(normData, "normalized_data.csv", row.names = TRUE)
```

> **Why gene-length normalization?** RNA-Seq read counts are proportional to both expression level AND gene length. Longer genes accumulate more reads simply due to size. Gene-length normalization removes this systematic bias before comparing expression across genes.

### Quantile Filtering

```r
# Remove bottom 25% lowest-expressed genes
fildata <- TCGAanalyze_Filtering(
  tabDF   = normData,
  method  = "quantile",
  qnt.cut = 0.25
)

dim(fildata)   # 34,539 genes × 494 samples
```

> **Why filter?** Low-expression genes are dominated by technical noise. Retaining them inflates the multiple testing burden and reduces statistical power. Removing the bottom quartile leaves 34,539 informative genes.

---

## 6. Differential Expression Analysis

### Separate Mutant and WT Expression Matrices

```r
# Barcodes by group
selectedBarcodesMutant <- subset(Metadata, IDH == "WT")$Barcode      # WT group
selectedBarcodesWT     <- subset(Metadata, IDH == "Mutant")$Barcode  # Mutant group

# Subset filtered matrix
mutantMatrix <- fildata[, selectedBarcodesMutant]
wtMatrix     <- fildata[, selectedBarcodesWT]

# Report sample counts
cat("WT samples:     ", ncol(mutantMatrix), "\n")
cat("Mutant samples: ", ncol(wtMatrix), "\n")
```

### Run edgeR-Based DEA via TCGAbiolinks

```r
dge_results <- TCGAanalyze_DEA(
  mat1       = mutantMatrix,
  mat2       = wtMatrix,
  Cond1type  = "Mutant",
  Cond2type  = "WT",
  pipeline   = "edgeR",    # Negative binomial GLM
  fdr.cut    = 0.01,       # Benjamini-Hochberg FDR threshold
  logFC.cut  = 1           # |log2FC| > 1 (2-fold change)
)

# Add expression level context for heatmap
results.level <- TCGAanalyze_LevelTab(
  dge_results,
  "Mutant", "WT",
  fildata[, selectedBarcodesMutant],
  fildata[, selectedBarcodesWT]
)
```

**DEA Parameters Explained:**

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `pipeline` | edgeR | Negative binomial model appropriate for count overdispersion |
| `fdr.cut` | 0.01 | Strict FDR to minimize false positives across 34K tests |
| `logFC.cut` | 1 | Minimum 2-fold change for biological relevance |

**Results:**

| Category | Count |
|----------|-------|
| Total DEGs | **5,916** |
| Upregulated in IDH-Mutant | **1,681** |
| Downregulated in IDH-Mutant | **4,235** |

---

## 7. Visualization

### Volcano Plot

```r
# Classification thresholds
logFC_cutoff  <- 1
pvalue_cutoff <- 0.01

# Assign significance categories
dge_results$threshold <- ifelse(
  dge_results$logFC > logFC_cutoff & dge_results$PValue < pvalue_cutoff,
  "Upregulated",
  ifelse(
    dge_results$logFC < -logFC_cutoff & dge_results$PValue < pvalue_cutoff,
    "Downregulated",
    "Not Significant"
  )
)

# Color coding
dge_results$color <- ifelse(dge_results$threshold == "Upregulated",   "red",
                     ifelse(dge_results$threshold == "Downregulated", "blue",
                                                                       "gray"))

# -log10 transform p-values
dge_results$log10_pvalue <- -log10(dge_results$PValue)

# Plot
plot(
  dge_results$logFC,
  dge_results$log10_pvalue,
  pch  = 16,
  col  = dge_results$color,
  xlab = "Log Fold Change (logFC)",
  ylab = "-log10 Adjusted P-value",
  main = "Volcano Plot: IDH Mutant vs Wildtype (LGG)",
  cex  = 1.2
)
abline(v = c(-1, 1),   lty = 2, col = "black")
abline(h = -log10(0.01), lty = 2, col = "black")

legend("topleft",
       legend = c("Upregulated", "Downregulated", "Not Significant"),
       col    = c("red", "blue", "gray"),
       pch    = 16)
```





<img width="1190" height="735" alt="Volcano_plot" src="https://github.com/user-attachments/assets/e33cd22b-ef12-48fb-ad94-bc3903462d60" />



**Figure 1:** Volcano plot of TCGA-LGG differential expression. Red points: upregulated in IDH-mutant (logFC > 1, p < 0.01). Blue points: downregulated in IDH-mutant. Dashed lines mark the thresholds.

---

### Heatmap

```r
# Extract expression for DEGs only
heat.data <- fildata[rownames(results.level), ]

# Sample group colors
numWT     <- length(subset(Metadata, IDH == "WT")$Barcode)
numMutant <- length(subset(Metadata, IDH == "Mutant")$Barcode)
cat("WT samples:    ", numWT,     "\n")
cat("Mutant samples:", numMutant, "\n")

# Color vector for sample annotation
ccols <- c(
  rep("red",  numWT),       # WT = red
  rep("blue", numMutant)    # Mutant = blue
)

# Generate heatmap with hierarchical clustering
heatmap.2(
  as.matrix(heat.data),
  col           = hcl.colors(10, "Blue-Red"),
  ColSideColors = ccols,
  scale         = "row",
  trace         = "none",
  labRow        = FALSE,
  labCol        = FALSE,
  key           = TRUE,
  keysize       = 1.5,
  density.info  = "none",
  main          = "DEG Expression Heatmap: IDH Mutant vs WT"
)
```



<img width="1000" height="800" alt="heatmap_output" src="https://github.com/user-attachments/assets/3410b4ab-2347-47b7-b982-fd73c8b07a51" />



**Figure 2:** Hierarchical clustering heatmap of 5,916 DEGs. Red column bar = WT samples; Blue column bar = Mutant samples. Row-scaled expression shows distinct separation between the two IDH status groups.

---

## 8. Functional Enrichment Analysis

### Separate Upregulated and Downregulated DEGs

```r
# Upregulated: positive logFC, significant FDR
upregulated_genes <- rownames(dge_results[
  dge_results$logFC > 1 & dge_results$FDR < 0.01, ])

# Downregulated: negative logFC, significant FDR
downregulated_genes <- rownames(dge_results[
  dge_results$logFC < -1 & dge_results$FDR < 0.01, ])

cat("Upregulated genes:  ", length(upregulated_genes), "\n")
cat("Downregulated genes:", length(downregulated_genes), "\n")
```

### GO Enrichment via TCGAbiolinks

```r
# Enrichment analysis — Upregulated genes (Biological Process)
ansUp <- TCGAanalyze_EAcomplete(
  TFname     = "upregulated",
  RegulonList = upregulated_genes
)

# Enrichment analysis — Downregulated genes
ansDown <- TCGAanalyze_EAcomplete(
  TFname      = "downregulated",
  RegulonList = downregulated_genes
)

# Barplot of top enriched GO:BP terms — Upregulated
TCGAvisualize_EAbarplot(
  tf          = rownames(ansUp$ResBP),
  GOBPTab     = ansUp$ResBP,
  nBar        = 10,
  title       = "GO Biological Process — Upregulated Genes",
  xlab        = "-log10(FDR)",
  ylab        = "GO Terms",
  filename    = "upregulated_GO_BP.pdf"
)


<img width="9000" height="4500" alt="Upregulated_EA" src="https://github.com/user-attachments/assets/1a110888-361b-434e-ab9f-faa7172bc5a1" />

Figure 3: Functional enrichment analysis of upregulated genes for low-grade gliomas (LGG).

---


# # Barplot of top enriched GO:BP terms — Downregulated
TCGAvisualize_EAbarplot(
  tf          = rownames(ansDown$ResBP),
  GOBPTab     = ansDown$ResBP,
  nBar        = 10,
  title       = "GO Biological Process — Downregulated Genes",
  xlab        = "-log10(FDR)",
  ylab        = "GO Terms",
  filename    = "downregulated_GO_BP.pdf"
)


<img width="9000" height="4500" alt="Downregulated_EA" src="https://github.com/user-attachments/assets/a747bb89-b270-4b0e-8fe5-32508cde6f2a" />


Figure 4: Functional enrichment analysis of downregulated genes for low-grade gliomas (LGG).

```

---

## 9. Results Summary

### Differentially Expressed Genes

| Metric | Value |
|--------|-------|
| Input genes (post-filter) | 34,539 |
| Significant DEGs (FDR < 0.01, \|logFC\| > 1) | **5,916** |
| Upregulated in IDH-Mutant | **1,681** |
| Downregulated in IDH-Mutant | **4,235** |
| Ratio (Down:Up) | ~2.5:1 |

### Top Enriched Pathways — Upregulated Genes

| GO Term | Biological Significance |
|---------|------------------------|
| Synaptic transmission | Altered neuronal signaling in IDH-mutant gliomas |
| Cell-cell signaling | Modified tumor-microenvironment communication |
| Homophilic cell adhesion | Potential changes in tumor cohesion |
| Transmission of nerve impulse | Neurophysiological rewiring |

### Top Enriched Pathways — Downregulated Genes

| GO Term | Biological Significance |
|---------|------------------------|
| Anterior/posterior pattern formation | Epigenetic silencing of developmental regulators |
| Skeletal system development | Downregulation of mesenchymal differentiation genes |
| Embryonic morphogenesis | IDH-mutant hypermethylation silencing developmental genes |
| Regionalization | Suppression of CNS patterning programs |

---

## 10. Interpretation

### Biological Context

The IDH mutation introduces a neomorphic enzymatic activity that produces 2-hydroxyglutarate (2-HG). This oncometabolite inhibits α-ketoglutarate-dependent dioxygenases, including TET methylcytosine dioxygenases and histone demethylases. The result is **genome-wide CpG island hypermethylation** (the G-CIMP phenotype), which silences many developmental and morphogenetic gene programs — precisely reflected in the downregulated GO terms observed here.

The **upregulated neural signaling** pathways may reflect the tumor's origin from glial progenitor cells retaining or acquiring neuronal communication properties, or represent compensatory signaling in response to IDH-mediated metabolic reprogramming.

### Why IDH-Mutant LGG Has Better Prognosis

The enrichment for developmental pathway silencing (downregulated genes) combined with preserved neural signaling (upregulated) is consistent with the well-established better clinical outcome of IDH-mutant LGG compared to IDH-wildtype glioblastoma — which tends to have more aggressive, mesenchymal, and invasive gene expression patterns.

### Analytical Strengths & Limitations

| Aspect | Detail |
|--------|--------|
| ✅ Strength | Gene-length normalization corrects length bias before DE |
| ✅ Strength | Strict FDR threshold (0.01) reduces false positives |
| ✅ Strength | edgeR negative binomial model appropriate for count data |
| ⚠️ Limitation | Class imbalance: 460 Mutant vs 34 WT |
| ⚠️ Limitation | Single-cohort analysis (TCGA only) |
| ⚠️ Limitation | Bulk RNA-Seq averages over cellular heterogeneity |

---

## Software Versions

| Package | Version | Role |
|---------|---------|------|
| R | ≥ 4.0 | Analysis environment |
| TCGAbiolinks | 2.x | Data access, normalization, DE |
| edgeR | 3.x | Differential expression GLM |
| EDASeq | 2.x | Normalization support |
| biomaRt | 2.x | Gene annotation |
| gplots | 3.x | Heatmap visualization |
| dplyr | 1.x | Data manipulation |
