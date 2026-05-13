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

<p align="center">
<svg width="680" height="860" viewBox="0 0 680 860" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <marker id="arr" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M2 1L8 5L2 9" fill="none" stroke="#888780" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arr-blue" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M2 1L8 5L2 9" fill="none" stroke="#185FA5" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arr-purple" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M2 1L8 5L2 9" fill="none" stroke="#534AB7" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arr-coral" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M2 1L8 5L2 9" fill="none" stroke="#993C1D" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
    <marker id="arr-amber" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M2 1L8 5L2 9" fill="none" stroke="#854F0B" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
  </defs>

  <!-- Step 01: Data source (teal) -->
  <rect x="150" y="20" width="380" height="56" rx="10" fill="#E1F5EE" stroke="#0F6E56" stroke-width="1"/>
  <text x="340" y="43" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#085041">TCGA-LUAD Dataset</text>
  <text x="340" y="62" text-anchor="middle" font-family="monospace" font-size="11" fill="#0F6E56">RNA-Seq · STAR-Counts · GDC Portal</text>
  <text x="134" y="52" text-anchor="middle" font-family="monospace" font-size="10" fill="#0F6E56">01</text>

  <line x1="340" y1="76" x2="340" y2="104" stroke="#0F6E56" stroke-width="1.5" marker-end="url(#arr)"/>

  <!-- Step 02: Query & Download (blue) -->
  <rect x="150" y="104" width="380" height="56" rx="10" fill="#E6F1FB" stroke="#185FA5" stroke-width="1"/>
  <text x="340" y="127" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#0C447C">Query &amp; Download</text>
  <text x="340" y="146" text-anchor="middle" font-family="monospace" font-size="11" fill="#185FA5">GDCquery → GDCdownload → GDCprepare</text>
  <text x="134" y="136" text-anchor="middle" font-family="monospace" font-size="10" fill="#185FA5">02</text>

  <line x1="340" y1="160" x2="340" y2="188" stroke="#185FA5" stroke-width="1.5" marker-end="url(#arr-blue)"/>

  <!-- Step 03: Metadata (blue) -->
  <rect x="150" y="188" width="380" height="56" rx="10" fill="#E6F1FB" stroke="#185FA5" stroke-width="1"/>
  <text x="340" y="211" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#0C447C">Create Metadata</text>
  <text x="340" y="230" text-anchor="middle" font-family="monospace" font-size="11" fill="#185FA5">Barcode · Race · Sample type · Sample ID</text>
  <text x="134" y="220" text-anchor="middle" font-family="monospace" font-size="10" fill="#185FA5">03</text>

  <line x1="340" y1="244" x2="340" y2="272" stroke="#185FA5" stroke-width="1.5" marker-end="url(#arr-blue)"/>

  <!-- Step 04: Downsize (purple) -->
  <rect x="150" y="272" width="380" height="56" rx="10" fill="#EEEDFE" stroke="#534AB7" stroke-width="1"/>
  <text x="340" y="295" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#3C3489">Downsize Dataset</text>
  <text x="340" y="314" text-anchor="middle" font-family="monospace" font-size="11" fill="#534AB7">20 Primary Tumor + 20 Solid Tissue Normal</text>
  <text x="134" y="304" text-anchor="middle" font-family="monospace" font-size="10" fill="#534AB7">04</text>

  <!-- Split to two purple boxes -->
  <path d="M340 328 L340 342 L200 342 L200 356" fill="none" stroke="#534AB7" stroke-width="1.5" marker-end="url(#arr-purple)"/>
  <path d="M340 328 L340 342 L480 342 L480 356" fill="none" stroke="#534AB7" stroke-width="1.5" marker-end="url(#arr-purple)"/>

  <!-- Step 05a: Normalization (purple) -->
  <rect x="80" y="356" width="236" height="52" rx="10" fill="#EEEDFE" stroke="#534AB7" stroke-width="1"/>
  <text x="198" y="378" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#3C3489">Normalization</text>
  <text x="198" y="396" text-anchor="middle" font-family="monospace" font-size="11" fill="#534AB7">geneLength method</text>
  <text x="64" y="385" text-anchor="middle" font-family="monospace" font-size="10" fill="#534AB7">05</text>

  <!-- Step 05b: Filtering (purple) -->
  <rect x="364" y="356" width="236" height="52" rx="10" fill="#EEEDFE" stroke="#534AB7" stroke-width="1"/>
  <text x="482" y="378" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#3C3489">Filtering</text>
  <text x="482" y="396" text-anchor="middle" font-family="monospace" font-size="11" fill="#534AB7">Quantile cut = 0.25</text>

  <!-- Merge back -->
  <path d="M198 408 L198 442 L340 442 L340 464" fill="none" stroke="#534AB7" stroke-width="1.5" marker-end="url(#arr-purple)"/>
  <path d="M482 408 L482 442 L340 442" fill="none" stroke="#534AB7" stroke-width="1.5"/>

  <!-- Step 06: DEA (coral) -->
  <rect x="130" y="464" width="420" height="68" rx="10" fill="#FAECE7" stroke="#993C1D" stroke-width="1"/>
  <text x="340" y="487" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#712B13">Differential Expression Analysis</text>
  <text x="340" y="506" text-anchor="middle" font-family="monospace" font-size="11" fill="#993C1D">edgeR · FDR &lt; 0.01 · |logFC| &gt; 1</text>
  <text x="340" y="523" text-anchor="middle" font-family="monospace" font-size="11" fill="#993C1D">↑ 3,277 upregulated   ↓ 6,357 downregulated</text>
  <text x="114" y="501" text-anchor="middle" font-family="monospace" font-size="10" fill="#993C1D">06</text>

  <!-- Split to two coral boxes -->
  <path d="M340 532 L340 546 L198 546 L198 560" fill="none" stroke="#993C1D" stroke-width="1.5" marker-end="url(#arr-coral)"/>
  <path d="M340 532 L340 546 L482 546 L482 560" fill="none" stroke="#993C1D" stroke-width="1.5" marker-end="url(#arr-coral)"/>

  <!-- Step 07a: Volcano (coral) -->
  <rect x="80" y="560" width="236" height="52" rx="10" fill="#FAECE7" stroke="#993C1D" stroke-width="1"/>
  <text x="198" y="582" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#712B13">Volcano Plot</text>
  <text x="198" y="600" text-anchor="middle" font-family="monospace" font-size="11" fill="#993C1D">logFC vs –log10 p-value</text>
  <text x="64" y="589" text-anchor="middle" font-family="monospace" font-size="10" fill="#993C1D">07</text>

  <!-- Step 07b: Heatmap (coral) -->
  <rect x="364" y="560" width="236" height="52" rx="10" fill="#FAECE7" stroke="#993C1D" stroke-width="1"/>
  <text x="482" y="582" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#712B13">Heatmap</text>
  <text x="482" y="600" text-anchor="middle" font-family="monospace" font-size="11" fill="#993C1D">Column clustering · Blue-Red 2</text>

  <!-- Merge back -->
  <path d="M198 612 L198 646 L340 646 L340 668" fill="none" stroke="#993C1D" stroke-width="1.5" marker-end="url(#arr-coral)"/>
  <path d="M482 612 L482 646 L340 646" fill="none" stroke="#993C1D" stroke-width="1.5"/>

  <!-- Step 08: Gene Annotation (amber) -->
  <rect x="150" y="668" width="380" height="56" rx="10" fill="#FAEEDA" stroke="#854F0B" stroke-width="1"/>
  <text x="340" y="691" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#633806">Gene Annotation</text>
  <text x="340" y="710" text-anchor="middle" font-family="monospace" font-size="11" fill="#854F0B">BioMart · Ensembl ID → HGNC Symbol</text>
  <text x="134" y="700" text-anchor="middle" font-family="monospace" font-size="10" fill="#854F0B">08</text>

  <line x1="340" y1="724" x2="340" y2="752" stroke="#854F0B" stroke-width="1.5" marker-end="url(#arr-amber)"/>

  <!-- Step 09: Enrichment (amber) -->
  <rect x="150" y="752" width="380" height="68" rx="10" fill="#FAEEDA" stroke="#854F0B" stroke-width="1"/>
  <text x="340" y="775" text-anchor="middle" font-family="monospace" font-size="13" font-weight="bold" fill="#633806">Functional Enrichment Analysis</text>
  <text x="340" y="794" text-anchor="middle" font-family="monospace" font-size="11" fill="#854F0B">GO: BP · CC · MF + KEGG Pathways</text>
  <text x="340" y="811" text-anchor="middle" font-family="monospace" font-size="11" fill="#854F0B">TCGAbiolinks · Barplot visualization</text>
  <text x="134" y="788" text-anchor="middle" font-family="monospace" font-size="10" fill="#854F0B">09</text>

  <!-- Legend -->
  <rect x="150" y="836" width="12" height="12" rx="2" fill="#E1F5EE" stroke="#0F6E56" stroke-width="1"/>
  <text x="168" y="847" font-family="monospace" font-size="10" fill="#555">Data source</text>
  <rect x="240" y="836" width="12" height="12" rx="2" fill="#E6F1FB" stroke="#185FA5" stroke-width="1"/>
  <text x="258" y="847" font-family="monospace" font-size="10" fill="#555">Preprocessing</text>
  <rect x="340" y="836" width="12" height="12" rx="2" fill="#EEEDFE" stroke="#534AB7" stroke-width="1"/>
  <text x="358" y="847" font-family="monospace" font-size="10" fill="#555">Normalization</text>
  <rect x="440" y="836" width="12" height="12" rx="2" fill="#FAECE7" stroke="#993C1D" stroke-width="1"/>
  <text x="458" y="847" font-family="monospace" font-size="10" fill="#555">Analysis</text>
  <rect x="510" y="836" width="12" height="12" rx="2" fill="#FAEEDA" stroke="#854F0B" stroke-width="1"/>
  <text x="528" y="847" font-family="monospace" font-size="10" fill="#555">Annotation</text>
</svg>
</p>

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

<p align="center">
  <img src="./Images/Volcano_Plot.png" width="650" alt="Volcano Plot showing upregulated (red) and downregulated (blue) genes in LUAD"/>
</p>

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

<p align="center">
  <img src="./Images/Heatmap.png" width="700" alt="Heatmap showing clear separation between tumor and normal tissue samples"/>
</p>

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

<p align="center">
  <img src="./Images/Upregulated_genes_EA.png" width="700" alt="GO and KEGG enrichment barplots for upregulated genes"/>
</p>

> **Figure 3.** Functional enrichment analysis of 3,277 upregulated genes across four categories: GO Biological Process, GO Cellular Component, GO Molecular Function, and KEGG Pathways. Top enriched biological processes include **cell adhesion**, **biological adhesion**, and **behavior** — pathways directly associated with tumor invasion, migration, and metastatic spread in LUAD. Top KEGG pathways include agranulocyte/granulocyte adhesion and diapedesis and GPCR signaling, consistent with immune evasion and pro-tumorigenic signaling. Bar length represents gene ratio; color depth reflects −log₁₀(FDR).

---

### Figure 4 — Functional Enrichment of Downregulated Genes (n = 6,357)

<p align="center">
  <img src="./Images/Downregulated_genes_EA.png" width="700" alt="GO and KEGG enrichment barplots for downregulated genes"/>
</p>

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
