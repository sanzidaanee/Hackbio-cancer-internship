# Differential Expression Analysis — Technical Report
## Glioma IDH Status Classification from TCGA-LGG RNA-Seq Data

**HackBio Cancer Bioinformatics Internship — Stage 4**  
**Language:** R &nbsp;|&nbsp; **Pipeline:** TCGAbiolinks + edgeR

---

## Table of Contents

1. [DEA Pipeline Workflow](#1-dea-pipeline-workflow)
2. [Environment Setup](#2-environment-setup)
3. [Step-by-Step Code Walkthrough](#3-step-by-step-code-walkthrough)
   - [Step 1 — Data Acquisition](#step-1--data-acquisition)
   - [Step 2 — Metadata Construction & IDH Status](#step-2--metadata-construction--idh-status)
   - [Step 3 — Normalization & Filtering](#step-3--normalization--filtering)
   - [Step 4 — Differential Expression Analysis](#step-4--differential-expression-analysis)
   - [Step 5 — Volcano Plot](#step-5--volcano-plot)
   - [Step 6 — Heatmap](#step-6--heatmap)
   - [Step 7 — Functional Enrichment Analysis](#step-7--functional-enrichment-analysis)
4. [Results](#4-results)
5. [Interpretation](#5-interpretation)
6. [Software Versions](#6-software-versions)

---

## 1. DEA Pipeline Workflow

![DEA Pipeline Workflow](./dea_pipeline_workflow.png)

The pipeline begins with raw TCGA-LGG RNA-Seq data and flows through seven structured steps: data acquisition, metadata construction, normalization, differential expression analysis, volcano plot visualization, heatmap clustering, and GO functional enrichment. The workflow splits at enrichment to examine upregulated and downregulated gene sets separately.

---

## 2. Environment Setup

### Install Required Packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")        # TCGA data access & analysis
BiocManager::install("edgeR")               # Differential expression
BiocManager::install("EDASeq")              # Normalization support
BiocManager::install("SummarizedExperiment")# Data structures
BiocManager::install("biomaRt")             # Gene ID annotation

install.packages("dplyr")   # Data manipulation
install.packages("gplots")  # Heatmap generation
```

### Load Libraries

```r
library(TCGAbiolinks)         # Accessing and querying TCGA data
library(edgeR)                # Differential expression analysis
library(EDASeq)               # Exploratory data analysis and normalization
library(SummarizedExperiment) # Access assay data
library(biomaRt)              # Access BioMart databases
library(dplyr)                # Manage and manipulate data
library(gplots)               # Data visualization and heatmap generation
```

---

## 3. Step-by-Step Code Walkthrough

---

### Step 1 — Data Acquisition

Query and download the TCGA-LGG RNA-Seq dataset alongside IDH mutation data.

```r
# Browse available TCGA projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-LGG')

# ── Build query for LGG gene expression data ─────────────────────────────────
lggQ <- GDCquery(
  project               = 'TCGA-LGG',
  data.category         = 'Transcriptome Profiling',
  experimental.strategy = "RNA-Seq",
  workflow.type         = "STAR - Counts",
  access                = "open",
  data.type             = "Gene Expression Quantification"
)

# ── Download and prepare into SummarizedExperiment object ────────────────────
GDCdownload(lggQ)
lgg.data <- GDCprepare(lggQ)

# Inspect data structure
head(lgg.data)
colnames(lgg.data)
nrow(lgg.data)   # 60,660 genes
ncol(lgg.data)   # 534 samples
```

```r
# ── Download IDH mutation data ────────────────────────────────────────────────
mutation_query <- GDCquery(
  project       = "TCGA-LGG",
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(mutation_query)
mutation_data <- GDCprepare(mutation_query)

# Filter for IDH1 and IDH2 mutations only
idh_mutations <- mutation_data[mutation_data$Hugo_Symbol %in% c("IDH1", "IDH2"), ]
```

**Output:** SummarizedExperiment object with 60,660 genes across 534 LGG patient samples.

---

### Step 2 — Metadata Construction & IDH Status

Extract clinical variables and align IDH mutation labels with RNA-Seq barcodes.

```r
# ── Explore metadata columns available ───────────────────────────────────────
colData(lgg.data)
col_df <- as.data.frame(colData(lgg.data))
head(col_df)

# Inspect key fields
lgg.data$barcode
lgg.data$paper_IDH.status           # Mutant / WT / NA
lgg.data$paper_TERT.expression.status
lgg.data$paper_ATRX.status

# ── Build metadata dataframe ──────────────────────────────────────────────────
Metadata <- data.frame(
  "Barcode"               = lgg.data$barcode,
  "Tumor_Sample_Barcode"  = lgg.data$bcr_patient_barcode,
  "tumor_type"            = lgg.data$tumor_descriptor,
  "sample_id"             = lgg.data$sample_id,
  "sample"                = lgg.data$sample_type,
  "gender"                = lgg.data$gender,
  "IDH"                   = lgg.data$paper_IDH.status,        # KEY column
  "Mutation"              = lgg.data$paper_Mutation.Count,
  "TERT"                  = lgg.data$paper_TERT.expression.status,
  "TERT_status"           = lgg.data$paper_TERT.promoter.status,
  "ATRX_status"           = lgg.data$paper_ATRX.status,
  "RNA_cluster"           = lgg.data$paper_IDH.specific.RNA.Expression.Cluster,
  "Random_Forest_cluster" = lgg.data$paper_Random.Forest.Sturm.Cluster
)

# Save metadata
write.csv(Metadata, "TCGA_LGG_metadata.csv", row.names = FALSE)
```

```r
# ── Align barcodes between RNA-Seq and mutation data ─────────────────────────
# Standardise to first 12 characters (patient-level ID)
idh_mutations$Tumor_Sample_Barcode <- substr(idh_mutations$Tumor_Sample_Barcode, 1, 12)
Metadata$Tumor_Sample_Barcode      <- substr(Metadata$Tumor_Sample_Barcode,      1, 12)

# Merge datasets on patient ID
merged_data <- merge(Metadata, idh_mutations, by = "Tumor_Sample_Barcode")

# Separate by IDH gene
idh1_data <- merged_data[merged_data$Hugo_Symbol == "IDH1", ]
idh2_data <- merged_data[merged_data$Hugo_Symbol == "IDH2", ]

# Separate by IDH mutation status
mutant_samples   <- merged_data[merged_data$IDH == "Mutant", ]
wildtype_samples <- merged_data[merged_data$IDH == "WT",     ]
```

**Output:** `TCGA_LGG_metadata.csv` — 534 rows with IDH status, TERT, ATRX, and cluster labels per sample.

---

### Step 3 — Normalization & Filtering

Extract raw counts, apply gene-length normalization, and remove low-expression genes.

```r
# ── Extract raw count matrix (unstranded) ─────────────────────────────────────
lgg.raw.data <- assays(lgg.data)
dim(lgg.raw.data$unstranded)   # 60,660 × 534

# ── Select only samples with a known IDH status ───────────────────────────────
selectedBarcodes <- c(
  subset(Metadata, IDH == "WT")$Barcode,
  subset(Metadata, IDH == "Mutant")$Barcode
)

selectedData <- lgg.raw.data$unstranded[, selectedBarcodes]
dim(selectedData)   # 60,660 × 494 (34 WT + 460 Mutant)

# Validate data
if (is.null(selectedData) || nrow(selectedData) == 0 || ncol(selectedData) == 0) {
  stop("Error: selectedData is empty. Check barcode alignment.")
} else {
  print(paste("Data OK — dimensions:", paste(dim(selectedData), collapse=" x ")))
}
```

```r
# ── Gene-length normalization ─────────────────────────────────────────────────
# Load gene length reference bundled with TCGAbiolinks
data(geneInfoHT, package = "TCGAbiolinks")

if (is.null(geneInfoHT) || nrow(geneInfoHT) == 0) {
  stop("geneInfoHT failed to load.")
} else {
  print("geneInfoHT loaded OK")
  head(geneInfoHT)
}

normData <- TCGAanalyze_Normalization(
  tabDF    = selectedData,
  geneInfo = geneInfoHT,
  method   = "geneLength"    # Corrects read counts for transcript length
)

# Save normalized matrix for downstream ML use
write.csv(normData, "normalized_data.csv", row.names = TRUE)
```

```r
# ── Quantile filtering: remove bottom 25% lowest-expressed genes ──────────────
fildata <- TCGAanalyze_Filtering(
  tabDF   = normData,
  method  = "quantile",
  qnt.cut = 0.25
)

dim(fildata)   # 34,539 genes × 494 samples
```

> **Why gene-length normalization?** RNA-Seq counts increase with gene length — a 10 kb gene accumulates more reads than a 1 kb gene even at equal expression. `geneLength` normalization removes this systematic bias before statistical testing.

> **Why quantile filtering?** Genes in the lowest 25% of expression are dominated by technical noise. Removing them reduces multiple-testing burden from 60,660 to 34,539 genes without losing biologically relevant signal.

---

### Step 4 — Differential Expression Analysis

Fit an edgeR negative binomial model to test IDH Mutant vs Wild-Type.

```r
# ── Assign barcodes to each condition ─────────────────────────────────────────
selectedBarcodesMutant <- subset(Metadata, IDH == "WT")$Barcode      # WT group
selectedBarcodesWT     <- subset(Metadata, IDH == "Mutant")$Barcode  # Mutant group

# Subset filtered expression matrix
mutantMatrix <- fildata[, selectedBarcodesMutant]
wtMatrix     <- fildata[, selectedBarcodesWT]

cat("WT samples:     ", ncol(mutantMatrix), "\n")
cat("Mutant samples: ", ncol(wtMatrix),     "\n")
```

```r
# ── Run differential expression analysis (edgeR) ──────────────────────────────
dge_results <- TCGAanalyze_DEA(
  mat1       = mutantMatrix,
  mat2       = wtMatrix,
  Cond1type  = "Mutant",
  Cond2type  = "WT",
  pipeline   = "edgeR",    # Negative binomial GLM
  fdr.cut    = 0.01,       # Benjamini-Hochberg FDR threshold
  logFC.cut  = 1           # Minimum |log2FC| = 1  (2-fold change)
)

# Add expression level context for heatmap
results.level <- TCGAanalyze_LevelTab(
  dge_results, "Mutant", "WT",
  fildata[, selectedBarcodesMutant],
  fildata[, selectedBarcodesWT]
)

head(results.level)
dim(results.level)   # 5,916 significant DEGs
```

**DEA Parameters:**

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `pipeline` | edgeR | Negative binomial GLM for overdispersed count data |
| `fdr.cut` | 0.01 | Strict FDR to minimize false discoveries across 34K tests |
| `logFC.cut` | 1 | Requires minimum 2-fold expression difference |
| Multiple testing | Benjamini-Hochberg | Standard FDR correction for genomics |

---

### Step 5 — Volcano Plot

Visualize the global landscape of differential expression.

```r
# ── Set significance thresholds ───────────────────────────────────────────────
logFC_cutoff  <- 1     # log2 fold change threshold
pvalue_cutoff <- 0.01  # adjusted p-value threshold

# ── Classify each gene ────────────────────────────────────────────────────────
dge_results$threshold <- ifelse(
  dge_results$logFC > logFC_cutoff  & dge_results$PValue < pvalue_cutoff, "Upregulated",
  ifelse(
    dge_results$logFC < -logFC_cutoff & dge_results$PValue < pvalue_cutoff, "Downregulated",
    "Not Significant"
  )
)

# Assign plot colors
dge_results$color <- ifelse(dge_results$threshold == "Upregulated",   "red",
                     ifelse(dge_results$threshold == "Downregulated", "blue", "gray"))

# Transform p-values for y-axis
dge_results$log10_pvalue <- -log10(dge_results$PValue)

# ── Generate volcano plot ──────────────────────────────────────────────────────
plot(
  dge_results$logFC,
  dge_results$log10_pvalue,
  pch  = 16,
  col  = dge_results$color,
  xlab = "Log Fold Change (logFC)",
  ylab = "-log10 Adjusted P-value",
  main = "Volcano Plot: Upregulated vs Downregulated",
  cex  = 1.2
)
abline(v = c(-logFC_cutoff, logFC_cutoff), lty = 2, col = "black")
abline(h = -log10(pvalue_cutoff),          lty = 2, col = "black")

legend("topleft",
       legend = c("Upregulated", "Downregulated", "Not Significant"),
       col    = c("red", "blue", "gray"),
       pch    = 16)
```

---

### Step 6 — Heatmap

Generate a hierarchical clustering heatmap of all 5,916 DEGs.

```r
# ── Prepare expression matrix for DEGs ────────────────────────────────────────
heat.data <- fildata[rownames(results.level), ]

dim(heat.data)   # 5,916 genes × 494 samples
ncol(heat.data)

# ── Count samples per group ───────────────────────────────────────────────────
numWT     <- length(subset(Metadata, IDH == "WT")$Barcode)
numMutant <- length(subset(Metadata, IDH == "Mutant")$Barcode)
cat("WT samples:    ", numWT,     "\n")   # 94
cat("Mutant samples:", numMutant, "\n")   # 419

# ── Assign column side colors ─────────────────────────────────────────────────
mutation.type <- c(rep("Mutant", 419), rep("WT", 94))
ccodes        <- c()
for (i in mutation.type) {
  if (i == "Mutant") {
    ccodes <- c(ccodes, "red")   # Red  = Mutant
  } else {
    ccodes <- c(ccodes, "blue")  # Blue = WT
  }
}

# Verify length matches columns
length(ccodes)   # Must equal ncol(heat.data)
```

```r
# ── Generate heatmap ──────────────────────────────────────────────────────────
png("heatmap_output.png", width = 1000, height = 800)

heatmap.2(
  x             = as.matrix(heat.data),
  col           = hcl.colors(10, palette = "Blue-Red 2"),
  Rowv          = TRUE,         # Enable row clustering
  Colv          = TRUE,         # Enable column clustering
  scale         = "row",        # Z-score scale per gene
  sepcolor      = "block",
  trace         = "none",
  key           = TRUE,
  dendrogram    = "both",       # Show dendrograms for rows and columns
  cexRow        = 0.5,
  cexCol        = 1,
  main          = "Heatmap of Mutant vs WT",
  na.color      = "black",
  ColSideColors = ccodes        # Red = Mutant, Blue = WT
)

dev.off()
```

---

### Step 7 — Functional Enrichment Analysis

Convert Ensembl IDs to gene symbols and run GO:BP enrichment on up- and downregulated gene sets.

```r
# ── Extract upregulated and downregulated gene sets ───────────────────────────
upreg_genes   <- rownames(subset(results.level, logFC >  1))
downreg_genes <- rownames(subset(results.level, logFC < -1))

cat("Upregulated genes:  ", length(upreg_genes),   "\n")   # 1,681
cat("Downregulated genes:", length(downreg_genes), "\n")   # 4,235
```

```r
# ── Convert Ensembl IDs to HGNC gene symbols via biomaRt ─────────────────────
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg_genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = upreg_genes,
  mart       = mart
)$hgnc_symbol

downreg_genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = downreg_genes,
  mart       = mart
)$hgnc_symbol

# Save gene lists
write.csv(upreg_genes,   "upregulated_genes.csv",   row.names = FALSE)
write.csv(downreg_genes, "downregulated_genes.csv", row.names = FALSE)
```

```r
# ── GO Biological Process enrichment ─────────────────────────────────────────
upreg_EA   <- TCGAanalyze_EAcomplete(TFname = "upregulated",   upreg_genes)
downreg_EA <- TCGAanalyze_EAcomplete(TFname = "downregulated", downreg_genes)
```

```r
# ── Visualize enrichment — Upregulated genes ──────────────────────────────────
TCGAvisualize_EAbarplot(
  tf       = rownames(upreg_EA$ResBP),
  GOBPTab  = upreg_EA$ResBP,
  GOCCTab  = upreg_EA$ResCC,
  GOMFTab  = upreg_EA$ResMF,
  PathTab  = upreg_EA$ResPat,
  nRGTab   = upreg_genes,
  nBar     = 5,
  text.size= 2,
  fig.width= 30,
  fig.height=15
)

# ── Visualize enrichment — Downregulated genes ────────────────────────────────
TCGAvisualize_EAbarplot(
  tf       = rownames(downreg_EA$ResBP),
  GOBPTab  = downreg_EA$ResBP,
  GOCCTab  = downreg_EA$ResCC,
  GOMFTab  = downreg_EA$ResMF,
  PathTab  = downreg_EA$ResPat,
  nRGTab   = downreg_genes,
  nBar     = 5,
  text.size= 2,
  fig.width= 30,
  fig.height=15
)
```

---

## 4. Results

### DEG Summary

| Category | Count |
|----------|-------|
| Input genes (raw) | 60,660 |
| Genes after filtering | 34,539 |
| **Total significant DEGs** | **5,916** |
| **Upregulated in IDH-Mutant** | **1,681** |
| **Downregulated in IDH-Mutant** | **4,235** |

---

### Figure 1 — Volcano Plot

![Volcano Plot: IDH Mutant vs Wild-Type LGG](./figures/Volcano_plot.png)

**Figure 1. Volcano plot of differentially expressed genes in TCGA-LGG (IDH Mutant vs Wild-Type).**  
The X-axis shows the log₂ fold change and the Y-axis shows the −log₁₀ adjusted p-value. **Red points** indicate significantly upregulated genes (logFC > 1, p < 0.01); **blue points** indicate significantly downregulated genes (logFC < −1, p < 0.01); **grey points** are not significant. Dashed vertical lines mark the ±1 logFC threshold and the dashed horizontal line marks FDR = 0.01. A total of 5,916 genes passed both thresholds.

---

### Figure 2 — Hierarchical Clustering Heatmap

![Heatmap of DEGs: IDH Mutant vs Wild-Type](./figures/heatmap_output.png)

**Figure 2. Hierarchical clustering heatmap of 5,916 differentially expressed genes across 494 TCGA-LGG samples.**  
Rows represent genes; columns represent patient samples. The **red column bar** indicates IDH-Mutant samples (n = 419); the **blue column bar** indicates Wild-Type samples (n = 94). Expression values are row-scaled (Z-score). Hierarchical clustering was applied to both rows and columns using the `heatmap.2()` function with the Blue-Red 2 colour palette. The heatmap reveals clear transcriptomic separation between IDH-Mutant and Wild-Type tumours, consistent with genome-wide G-CIMP hypermethylation in IDH-Mutant LGG.

---

### Figure 3 — GO Enrichment: Upregulated Genes

![GO Enrichment Barplot — Upregulated Genes](./figures/Upregulated_EA.png)

**Figure 3. Gene Ontology (Biological Process) enrichment barplot for 1,681 upregulated genes in IDH-Mutant LGG.**  
The top enriched GO:BP terms are shown ranked by −log₁₀(FDR). Upregulated pathways are predominantly associated with neural signalling and cell communication, including **synaptic transmission**, **cell-cell signalling**, **homophilic cell adhesion**, and **transmission of nerve impulse**. These findings suggest that IDH-Mutant tumours retain neuronal-like signalling properties, consistent with their glial progenitor cell origin and relatively better prognosis compared to IDH-Wildtype glioblastoma.

---

### Figure 4 — GO Enrichment: Downregulated Genes

![GO Enrichment Barplot — Downregulated Genes](./figures/Downregulated_genes_EA.png)

**Figure 4. Gene Ontology (Biological Process) enrichment barplot for 4,235 downregulated genes in IDH-Mutant LGG.**  
Downregulated pathways are dominated by developmental and morphogenetic processes, including **anterior/posterior pattern formation**, **skeletal system development**, **regionalization**, **pattern specification**, and **embryonic morphogenesis**. The silencing of these developmental gene programs is mechanistically explained by genome-wide CpG island hypermethylation (G-CIMP phenotype) driven by 2-hydroxyglutarate (2-HG), the oncometabolite produced by mutant IDH1/IDH2 enzymes.

---

## 5. Interpretation

### Biological Context

The IDH1/IDH2 mutations introduce a neomorphic enzyme activity that produces **2-hydroxyglutarate (2-HG)**. This oncometabolite competitively inhibits α-ketoglutarate-dependent dioxygenases — including TET methylcytosine dioxygenases and histone demethylases — resulting in:

- **Genome-wide CpG island hypermethylation** (the G-CIMP phenotype)
- **Silencing of developmental and morphogenetic gene programs** (reflected in downregulated GO terms)
- **Altered chromatin accessibility and transcription factor binding**

The upregulated neural signalling pathways likely reflect the tumour's origin from glial progenitor cells that retain neuronal communication properties, or represent compensatory signalling responses to IDH-mediated metabolic reprogramming.

### Why IDH-Mutant LGG Has Better Prognosis

The transcriptomic profile identified here — preserved neural communication (upregulated) combined with silenced invasive/mesenchymal developmental programs (downregulated) — is mechanistically consistent with the established clinical observation that IDH-Mutant LGG carries significantly better survival outcomes than IDH-Wildtype glioblastoma, which exhibits pro-invasive, mesenchymal gene expression signatures.

### Analytical Strengths & Limitations

| Aspect | Detail |
|--------|--------|
| ✅ Strength | Gene-length normalization removes transcript-size bias |
| ✅ Strength | Strict FDR = 0.01 minimises false positives across 34K genes |
| ✅ Strength | edgeR negative binomial model appropriate for count overdispersion |
| ⚠️ Limitation | Severe class imbalance: 460 Mutant vs 34 WT |
| ⚠️ Limitation | Single-cohort analysis — validation on CGGA or other datasets recommended |
| ⚠️ Limitation | Bulk RNA-Seq averages over tumour cellular heterogeneity |

---

## 6. Software Versions

| Package | Version | Role |
|---------|---------|------|
| R | ≥ 4.0 | Analysis environment |
| TCGAbiolinks | 2.x | Data access, normalization, DEA |
| edgeR | 3.x | Differential expression GLM |
| EDASeq | 2.x | Normalization support |
| biomaRt | 2.x | Ensembl ID to HGNC gene symbol conversion |
| gplots | 3.x | Heatmap visualization |
| dplyr | 1.x | Data manipulation |
| SummarizedExperiment | 1.x | Bioconductor data structure |
