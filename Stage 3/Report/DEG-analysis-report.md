# Differential Expression & Functional Enrichment Analysis

**HackBio Cancer Bioinformatics Internship — Stage 3**  
**Author:** Sanzida Akhter Anee  
**Dataset:** TCGA Lung Adenocarcinoma (LUAD)  
**Language:** R (RMarkdown)

---

## Overview

This pipeline performs differential gene expression (DEG) analysis between **Primary Tumor** and **Solid Tissue Normal** samples from the TCGA-LUAD dataset. It then annotates DEGs with HGNC gene symbols and runs functional enrichment analysis to reveal the biological pathways disrupted in lung adenocarcinoma.

---

## Pipeline Workflow

<img width="2351" height="3420" alt="deg_workflow_diagram" src="https://github.com/user-attachments/assets/0aa14e5d-a559-4f7b-a1ab-3a66d6ab9dc7" />



---

## Dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "TCGAbiolinks",
    "edgeR",
    "EDASeq",
    "SummarizedExperiment",
    "biomaRt"
))

install.packages(c("dplyr", "gplots"))
```

| Package | Purpose |
|---------|---------|
| `TCGAbiolinks` | Query, download, and prepare TCGA-LUAD data |
| `edgeR` | Statistical DEA using negative binomial model |
| `EDASeq` | Gene-length normalization of raw counts |
| `SummarizedExperiment` | Access assay slots from TCGA objects |
| `biomaRt` | Map Ensembl IDs to HGNC gene symbols |
| `dplyr` | Data manipulation and gene filtering |
| `gplots` | Heatmap generation via `heatmap.2` |

---

## Step-by-Step Code Walkthrough

### Step 1 — Query & Download TCGA-LUAD

```r
luadQ <- GDCquery(
    project               = 'TCGA-LUAD',
    data.category         = 'Transcriptome Profiling',
    experimental.strategy = "RNA-Seq",
    workflow.type         = "STAR - Counts",
    access                = "open",
    data.type             = "Gene Expression Quantification",
    sample.type           = c("Primary Tumor", "Solid Tissue Normal")
)
GDCdownload(luadQ)
luad.data <- GDCprepare(luadQ)
```

### Step 2 — Create Metadata & Downsize to 40 Samples

```r
Metadata <- data.frame(
    "barcode"    = luad.data$barcode,
    "race"       = luad.data$race,
    "tumor_type" = luad.data$tumor_descriptor,
    "sample"     = luad.data$sample_type,
    "sample_id"  = luad.data$sample_id
)

selectedBarcodes <- c(
    subset(Metadata, sample == "Primary Tumor")$barcode[1:20],
    subset(Metadata, sample == "Solid Tissue Normal")$barcode[1:20]
)
selectedData <- luad.raw.data$unstranded[, selectedBarcodes]
```

### Step 3 — Normalize & Filter

```r
normData <- TCGAanalyze_Normalization(
    tabDF = selectedData, geneInfo = geneInfoHT, method = "geneLength"
)
fildata <- TCGAanalyze_Filtering(
    tabDF = normData, method = "quantile", qnt.cut = 0.25
)
```

### Step 4 — Differential Expression Analysis

```r
results <- TCGAanalyze_DEA(
    mat1      = fildata[, selectedBarcodes[1:20]],
    mat2      = fildata[, selectedBarcodes[21:40]],
    Cond1type = "Primary Tumor",
    Cond2type = "Solid Tissue Normal",
    pipeline  = "edgeR",
    fdr.cut   = 0.01,
    logFC.cut = 1
)
```

| Threshold | Value | Result |
|-----------|-------|--------|
| FDR adjusted p-value | < 0.01 | — |
| Log₂ Fold Change | > 1 | **3,277 upregulated** |
| Log₂ Fold Change | < −1 | **6,357 downregulated** |

### Step 5 — Volcano Plot

```r
results$color <- ifelse(results$logFC >  1 & results$PValue < 0.01, "red",
                 ifelse(results$logFC < -1 & results$PValue < 0.01, "blue", "gray"))

plot(results$logFC, -log10(results$PValue),
     pch = 16, col = results$color,
     xlab = "Log Fold Change (logFC)",
     ylab = "-log10 Adjusted P-value",
     main = "Volcano Plot: Upregulated vs Downregulated")
```

---

## Results

### Figure 1 — Volcano Plot of Differential Expression (Primary Tumor vs Solid Tissue Normal)

<img width="1400" height="865" alt="Volcano Plot" src="https://github.com/user-attachments/assets/c89a3b32-4549-4ec1-a01f-37ed04945e21" />



> **Figure 1.** Volcano plot displaying differential gene expression between Primary Tumor (n = 20) and Solid Tissue Normal (n = 20) TCGA-LUAD samples. The x-axis shows log₂ fold change and the y-axis shows −log₁₀ adjusted p-value. **Red dots** (right, logFC > 1) represent 3,277 significantly upregulated genes in tumor tissue — candidate oncogenes and invasion drivers. **Blue dots** (left, logFC < −1) represent 6,357 downregulated genes — including potential tumor suppressors and cell cycle regulators. Significance threshold: FDR < 0.01, |logFC| > 1.

---

### Figure 2 — Heatmap of Significant DEGs Across All 40 Samples

```r
heatmap.2(
    x = as.matrix(heat.data),
    col = hcl.colors(10, palette = "Blue-Red 2"),
    Rowv = FALSE, Colv = TRUE, scale = "row",
    trace = "none", dendrogram = "col",
    ColSideColors = ccodes
)
```

<img width="1000" height="800" alt="Heatmap" src="https://github.com/user-attachments/assets/6cbb4c04-077b-46ed-bee4-9236509559a1" />


> **Figure 2.** Heatmap of all statistically significant DEGs (FDR < 0.01) across 40 samples. The dendrogram (top) clusters samples by expression similarity — **blue bars** denote Solid Tissue Normal; **red bars** denote Primary Tumor. The two groups separate into clearly distinct branches, confirming strong transcriptional differences between conditions. Red cells indicate high expression in that sample; blue cells indicate low expression. Genes expressed highly in tumor but lowly in normal are candidate oncogenes; the reverse pattern identifies candidate tumor suppressors.

---

### Step 6 — Gene Annotation via BioMart

```r
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg_genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = upreg_genes,
    mart       = mart
)$hgnc_symbol
```

### Step 7 — Functional Enrichment Analysis

```r
upreg_EA   <- TCGAanalyze_EAcomplete(TFname = "upregulated",   upreg_genes)
downreg_EA <- TCGAanalyze_EAcomplete(TFname = "downregulated", downreg_genes)

TCGAvisualize_EAbarplot(
    tf = rownames(upreg_EA$ResBP), GOBPTab = upreg_EA$ResBP,
    GOCCTab = upreg_EA$ResCC,      GOMFTab  = upreg_EA$ResMF,
    PathTab = upreg_EA$ResPat,     nRGTab   = upreg_genes,
    nBar = 5, fig.width = 30,      fig.height = 15
)
```

---

### Figure 3 — Functional Enrichment of Upregulated Genes (n = 3,277)


<img width="7500" height="3750" alt="Upregulated_genes_EA" src="https://github.com/user-attachments/assets/8c3c513d-d4b7-4f93-94b5-9b91186de0e6" />


> **Figure 3.** Functional enrichment analysis of 3,277 upregulated genes across four categories: GO Biological Process, GO Cellular Component, GO Molecular Function, and KEGG Pathways. Top enriched biological processes include **cell adhesion**, **biological adhesion**, and **behavior** — pathways directly associated with tumor invasion, migration, and metastatic spread in LUAD. Top KEGG pathways include agranulocyte/granulocyte adhesion and diapedesis and GPCR signaling, consistent with immune evasion and pro-tumorigenic signaling. Bar length represents gene ratio; color depth reflects −log₁₀(FDR).

---

### Figure 4 — Functional Enrichment of Downregulated Genes (n = 6,357)

<img width="7500" height="3750" alt="Downregulated_genes_EA" src="https://github.com/user-attachments/assets/085b810e-565a-4ae8-b8e7-c15ac08bc40c" />


> **Figure 4.** Functional enrichment analysis of 6,357 downregulated genes. Top enriched biological processes include **cell-cell signaling**, **nuclear division**, and **mitosis** — indicating disruption of normal cell cycle control and intercellular communication in LUAD. Top KEGG pathways include cell cycle control of chromosomal replication and mitotic roles of Polo-like kinase, consistent with loss of cell cycle checkpoints in cancer. Downregulation of these processes suggests dysregulated proliferation and impaired tissue homeostasis.

---

## Output Files

| File | Description |
|------|-------------|
| `TCGA_LUAD_metadata.csv` | Sample metadata (barcode, race, sample type) |
| `upregulated_genes.levels.csv` | Upregulated DEGs with fold-change values |
| `downregulated_genes.levels.csv` | Downregulated DEGs with fold-change values |
| `upregulated_genes_ei.csv` | Upregulated HGNC symbols after BioMart annotation |
| `downregulated_genes_ei.csv` | Downregulated HGNC symbols after BioMart annotation |
| `heatmap_output.png` | Heatmap of all significant DEGs |

---

## Links

| Resource | URL |
|----------|-----|
| GitHub Repository | https://github.com/sanzidaanee/Hackbio-cancer-internship/tree/main/Stage%203 |
| Full R Code | https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/Stage%203/Code/DEG.Rmd |
| Dataset | https://github.com/sanzidaanee/Hackbio-cancer-internship/tree/main/Stage%203/Data |

---

## Authors

Sanzida Akhter Anee (@Sanzida), Sk Arif (@arif_shaikh), Nada Ghozlan (@Nad1), Mennatallah Mohamed Ebrahim Mahmoud (@Mennatallah), Stéphie Raveloson (@StephieRav), Chidimma Nwaku (@Mma), Zaka Ullah (@Zaka), Idahosa Clinton (@doc_idahosa)
