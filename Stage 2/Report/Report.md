
# Gene Expression Data Visualization and Interpretation to Generate a Heatmap and Perform Downstream Functional Enrichment Analysis

Authors (@slack): Sanzida Akhter Anee (@Sanzida), Sk Arif (@arif_shaikh), Mina Zakaria (@mina_zakaria), Nada Ghozlan (@Nad1), Mennatallah Mohamed Ebrahim Mahmoud (@Mennatallah), Stéphie Raveloson (@StephieRav) 


## R Code link

https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/code/gene_expression_heatmap.Rmd


## Heatmaps for Gene Expression Analysis
### What is a Heatmap?


A heatmap is a common method of visualizing data. In the context of gene expression data, it displays the expression of many genes across many samples and helps to find patterns in gene expression data by aggregating data points and turning them into a simple visual representation making it easier to find a pattern [1]. 

In heatmaps, the data displayed in a grid, each row represents a gene and each column represents a sample. The changes in gene expression is represented by color and intensity of the boxes.   



### Gene Expression Dataset

Glioblastoma gene expression dataset with 500+ differentially expressed genes under different conditions.


Glioblastoma, also known as glioblastoma multiforme (GBM), is an aggressive type of cancer that occurs in the brain or spinal cord. The data represent the samples of cancer tissue from the patients with gene expression data for many different types of genes. 

### Link of Data

(https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/Stage%202/Data/glioblastoma.csv)

## Heatmap from Glioblastoma Dataset

In the heatmap, glioblastoma samples, with the choice of color palette (diverging vs. sequential) plays an important role in understanding the patterns of gene expression. 



## Part1: Data pre-processing
### Install and load Required Packages


```bash
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('biomaRt')
install.packages("gplots")
install.packages("dplyr")
install.packages("ggplot2")

```

#### Load required libraries

```bash
llibrary(biomaRt)
library(gplots)
library(dplyr)
library(ggplot2)
```


### Download data 
 ```bash
data <- read.csv ("/Users/sanzidaakhteranee/Documents/Internship_2024/HackBio/sTAGE 2/glioblastoma.csv", header = TRUE, row.names = 1 )
```


### Data structure  

 ```bash
head(data)
View(data)
colnames(data)       
nrow(data)             
ncol(data)
```



## Part 2: Heatmap Generation

 ### Generate heatmap   

```bash
  heatmap.2(as.matrix(data), trace = 'none')
```

![Fig 1](https://github.com/user-attachments/assets/4e2d14db-f108-44c8-931b-3c699b31af9a)








    

 
 ### Scaling

```bash
  heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE)
```

![Fig 2](https://github.com/user-attachments/assets/ed8d6cfa-11ed-4c3b-b4c5-9c45bf5d0a1f)








### Heatmap with diverging color palettes

```bash
  heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE,
          col=hcl.colors(100, palette = 'green-brown'))
```

![Fig 3](https://github.com/user-attachments/assets/249a73dd-6aac-4031-b746-3bcf49d7b7d3)




    


### Heatmap with sequential color palettes
```bash
  heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE,
          col=hcl.colors(100, palette = 'Blues3'))
```


![Fig 4](https://github.com/user-attachments/assets/15d33eb8-1136-4a70-a82c-6a1232a64056)





## Importance of Color Selection  in Gene Expression Heatmaps


Creating two color variants of the same heatmap by using diverging and sequential color palettes in gene expression data, with heatmap.2() function from the gplots package in R. The choice of color palette is important for interpreting gene expression data effectively to distinguish between different patterns in gene regulation. 

We generate heatmaps by using two color palettes:

- Diverging color palette (e.g., green to brown): represent both upregulation and downregulation in gene expression
- Sequential color palette (e.g., white to blue) : highlight the intensity of gene expression, irrespective of the direction of regulation



### Heatmap with Diverging Palette:

- This color palette is commonly used for distinguishing between two opposing conditions  to display both upregulated and downregulated genes, where one color represent higher expression (upregulation) and the opposite color represents lower expression (downregulation)

- This color palette (green to brown) can easily identify the contrasting gene expression patterns. Such as, overexpressed genes (positive fold change) can be represented in green, while underexpressed genes (negative fold change) are in brown, with white representing neutral expression

### Heatmap with Sequential Palette:

- This color palette (shades of blue), focuses on how strongly a gene is expressed with range of expression levels rather than it’s relative increase or decrease
- This palette is well-suited for identifying highly expressed genes without considering  up- or downregulated expression. It can be used when absolute expression values are more important than differential expression



## Heatmap Results Interpretation from Glioblastoma Dataset

- In diverging color palette, genes that are predominantly green in the glioblastoma columns are upregulated while  brown columns means these these genes are likely overexpressed in the glioblastoma condition
- Genes that are predominantly brown in the glioblastoma columns and green columns are likely underexpressed in glioblastoma condition
- This visual differentiation aids in identifying genes that may serve as potential biomarkers or pathways involved in glioblastoma progression or tumorigenesis
- In a sequential color palette, rows with dark blue values across multiple columns  indicate that the gene is highly expressed in these conditions while rows that are primarily white or light blue represent genes with low expression levels across samples
- The intensity of gene expression varies across conditions. Genes with darker shades (higher expression) in glioblastoma samples compared to lighter shades may represent genes that are overexpressed in glioblastoma condition. However, it doesn’t provide immediate information about whether the gene is upregulated or downregulated because it lacks a contrasting color for downregulation



## Part3: Heatmap Clustering

The heatmap can be combined with clustering methods which represent  group genes with or without samples together based on the similarity of their gene expression pattern. This is useful for identifying genes that are commonly regulated, or biological signatures associated with a particular disease condition

### Heatmap with Clustering of Genes (Rows Only)


```bash
  heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'row', 
          Colv = FALSE, Rowv = TRUE,
          col=hcl.colors(100, palette = 'green-brown'))
```

![rows only](https://github.com/user-attachments/assets/ed1786f1-4e1e-44b8-9830-d67ae6fb35e0)





### Heatmap with Clustering of Samples (Columns Only)

```bash
  heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE,
          col=hcl.colors(100, palette = 'green-brown'))
```

![Col only](https://github.com/user-attachments/assets/c578672c-e781-40a2-a39a-21a23a3bc19f)





### Heatmap with Clustering both Genes and Sample Together

```bash
  heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'both', 
          Colv = TRUE, Rowv = TRUE,
          col=hcl.colors(100, palette = 'green-brown'))

```

![both rows and col](https://github.com/user-attachments/assets/e5976f0b-d1e1-4954-be84-effafc5f28e6)





## Heatmap Clustering  Results Interpretation from Glioblastoma Dataset


- In a glioblastoma dataset, clusters of genes that are overexpressed in glioblastoma samples might be related to tumor growth or invasion, and their co-expression suggests they might be part of the same pathway
- Clustering of samples, all the glioblastoma samples group together in one cluster, suggesting that the glioblastoma samples have distinct gene expression profiles 
- Clustering both genes and samples together, shows  that a specific group of genes (upregulated in green) are consistently expressed across all glioblastoma samples that would be a  highlight of a potential set of glioblastoma-specific biomarkers
- Similarly, subgroups of glioblastoma samples cluster separately, indicating tumor heterogeneity, which could represent different genetic mutations or subtypes within glioblastoma.



## Part 4: Fold Change and P Value Calculation


### Selecting groups by index positions

```bash
group1<- c(1, 2,3,4,5)
group2<- c(6, 7,8,9,10)
```


### Group 1 & 2 from data 

```bash
group1_data <- data[, group1]
group2_data <- data[, group2]
```

### Calculate group mean

```bash
group1_mean <- rowMeans(group1_data)
group2_mean <- rowMeans(group2_data)
```

### Calculate fold change

```bash
fold_change <- (group2_mean)-(group1_mean)/group1_mean
logFC <- log2(fold_change)
```

### Calculate P value

```bash
  pvalues <- apply(data, 1, function(row) {
  t.test(row[1:5], row[6:10])$p.value
})
```


### Add fold change and p-values to the dataset

```bash
  results <- data.frame(Gene = rownames(data), logFC, pvalues)
```


## Part 5: Subset genes


### Set cutoff thresholds for fold change and p-value

```bash
logFC_up_cutoff <-  1   # Log2 fold change > 1 for upregulation
logFC_down_cutoff <- -1 # Log2 fold change < -1 for downregulation
pvalue_cutoff <- 0.05
```


### Subset significantly upregulated genes


```bash
upregulated_genes <- results %>%
  filter(logFC > logFC_up_cutoff & pvalues < pvalue_cutoff)

```


### Subset significantly downregulated genes

```bash
downregulated_genes <- results %>%
  filter(logFC < logFC_down_cutoff & pvalues < pvalue_cutoff)
```


### Print the results

```bash
print(upregulated_genes)
print(downregulated_genes)
```
    

### View the first few rows of the upregulated and downregulated genes

```bash
head(upregulated_genes)
head(downregulated_genes)
```


### Save upregulated genes to CSV

```bash
write.csv(upregulated_genes, "upregulated_genes.csv")
```


### Save downregulated genes to CSV

```bash
write.csv(downregulated_genes, "downregulated_genes.csv")
```

## Convert Ensembl IDs to gene id (upregulated genes) using biomaRt

BioMart is a powerful, web-based data management and analysis tool that integrates data from multiple biological databases, including Ensembl (genome data) and others


### Connect to the Ensembl BioMart database for human genes

```bash
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
```

### List of Ensembl IDs from upregulated genes data

```bash
ensembl_ids <- c("ENSG00000241945", "ENSG00000279104", "ENSG00000231107", 
                 "ENSG00000254092", "ENSG00000172236", "ENSG00000197253", 
                 "ENSG00000172116", "ENSG00000162598", "ENSG00000256193", 
                 "ENSG00000160183")
```
    

### Retrieve gene symbols for the Ensembl IDs

 ```bash
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                         filters = "ensembl_gene_id",
                         values = ensembl_ids,
                         mart = ensembl)
```

### View the results

```bash
print(gene_conversion)
```



## Functional Enrichment Analysis

 Functional enrichment analysis is a process to identify biological pathways, functions, or processes that are over-represented in a given set of genes that provide information into the biological mechanisms underlying a specific condition, disease progression or response to treatment.

We used the ShinyGO (Version: v0.741) tool with the GO biological process and P-value cutoff (FDR) is 0.5  to visualize the top 10 pathways. From here we choose to describe the top 3 enriched pathways based on  biological processes. 

![pathways](https://github.com/user-attachments/assets/7ac656e4-e38e-42d7-8b0e-c424fb60f2b9)








## Pathway Visualization

To visualize top 10 pathways for upregulated gene of glioblastoma dataset by using lollipop plot with scaling the points according to the negative log10 of the p-value.

<img width="930" alt="download" src="https://github.com/user-attachments/assets/ec0b3321-bf0a-4698-bc33-d795b0e9aa17">



## Pathway Visualization in R

```bash
#Input data from the table in the image
data <- data.frame(
  pathway = c("Ribosomal small subunit assembly", 
              "Maturation of SSU-rRNA from tricistronic rRNA", 
              "Proteolysis", 
              "Cellular sodium ion homeostasis", 
              "Maturation of SSU-rRNA", 
              "Ribosome assembly", 
              "Regulation of defense response to virus by virus", 
              "Sodium ion homeostasis", 
              "Ribosomal small subunit biogenesis", 
              "Extracellular matrix disassembly"),
  nGenes = c(1, 1, 3, 1, 1, 1, 1, 1, 1, 1),
  FDR = c(1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.7e-01, 1.7e-01),
  fold_enrichment = c(200, 88.4, 5.6, 172.7, 66.7, 59.4, 126.6, 65.5, 45.2, 40.4)
)

# Calculate -log10(FDR) for significance scaling
data$log_FDR <- -log10(data$FDR)

 

# Create the lollipop plot
ggplot(data, aes(y = reorder(pathway, nGenes), x = nGenes)) +  # Switch axes
  geom_segment(aes(y = reorder(pathway, nGenes), 
                   yend = reorder(pathway, nGenes), 
                   x = 0, 
                   xend = nGenes), color = "gray") +  # Lollipop stems
  geom_point(aes(size = log_FDR, color = fold_enrichment), alpha = 0.7) +  # Lollipop heads
  scale_size_continuous(range = c(3, 10)) +  # Adjust size based on -log10(FDR)
  scale_color_gradient(low = "lightblue", high = "darkblue") +  # Color based on fold enrichment
  labs(x = "Number of Genes", y = "Pathway", 
       size = "-log10(FDR)", color = "Fold Enrichment",
       title = "Pathways and Associated Genes with Significance") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  # Ensure y-axis labels are readable
```

![Rplot](https://github.com/user-attachments/assets/6b1d8145-a57d-48a2-a92b-4f0365ce318f)






## Top 3 Enriched Pathways

### 1. Ribosomal small subunit assembly


- Ribosomal small subunit assembly is an essential process within the broader pathway of ribosome biogenesis, which plays a crucial role in protein synthesis and the small ribosomal subunit, 40S in eukaryotes and 30S in prokaryotes, that is primarily responsible for decoding messenger RNA (mRNA) during translation[2, 3].
- In glioblastoma, an aggressive form of brain cancer, dysregulation of ribosome biogenesis and small subunit assembly is associated with enhanced cellular growth and proliferation, which is a hallmark of cancer,  and with the increase of ribosome biogenesis, the tumor cells gain the ability to produce proteins at an elevated rate, facilitating their rapid growth and survival in the tumor microenvironment [4].

### 2. Maturation of SSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)

- The Maturation of the small subunit ribosomal RNA (SSU-rRNA) from the tricistronic rRNA transcript is one of the vital biological process in ribosome biogenesis, specifically responsible for producing functional 18S rRNA, which is a component of the small ribosomal subunit and this process involves cleaving and modifying the precursor tricistronic transcript, which also encodes the 5.8S rRNA and large subunit rRNA (LSU-rRNA) [5,6]. 
- In cancer, such as glioblastoma, elevated ribosome biogenesis is frequently observed, where pathways involving rRNA processing and assembly are upregulated to sustain the high levels of protein synthesis needed for tumor growth and targeting these pathways, including the rRNA processing steps [7, 8], offers a potential therapeutic approach for inhibiting tumor progression [9].

### 3. Proteolysis

- Proteolysis involves  the breakdown of proteins into smaller polypeptides or amino acids through the action of enzymes called proteases which is essential for maintaining cellular homeostasis, regulating protein turnover, and enabling various signaling pathways [10]. 
- In cancer and other diseases, proteolysis plays critical roles in tumor progression, metastasis, and invasion by controlling the degradation of key regulatory proteins, extracellular matrix components, and other substrates [11].





# References

1. Wang, W., et al. (2013). Visualization of large-scale gene expression data: A survey. BMC Bioinformatics, 14(Suppl 8), S9. https://doi.org/10.1186/1471-2105-14-S8-S9.
2. Sloan, K. E., et al. (2017). The roles of snoRNAs in rRNA modification and ribosome assembly. Nature Reviews Molecular Cell Biology, 18(3), 190-205.
3. Peña, C., Hurt, E., & Panse, V. G. (2017). Eukaryotic ribosome assembly, transport and quality control. Nature Structural & Molecular Biology, 24(9), 689-699.
4. Morita, M., et al. (2021). Targeting ribosomal biogenesis in cancer. Nature Reviews Cancer, 21(3), 200-217.
5. Woolford, J. L., & Baserga, S. J. (2013). Ribosome biogenesis in the yeast Saccharomyces cerevisiae. Genetics, 195(3), 643-681.
6. Sloan, K. E., et al. (2017). The roles of snoRNAs in rRNA modification and ribosome assembly. Nature Reviews Molecular Cell Biology, 18(3), 190-205.

7. Hwang, S. P., & Denicourt, C. (2024). The impact of ribosome biogenesis in cancer: from proliferation to metastasis. NAR cancer, 6(2), zcae017.
8. Pelletier, J., Thomas, G., & Volarevic, S. (2018). Ribosome biogenesis in cancer: new players and therapeutic avenues. Nature Reviews Cancer, 18(1), 51-63. https://doi.org/10.1038/nrc.2017.104.
9. Zhang, S., et al. (2017). Targeting ribosome biogenesis in glioblastoma cells. Oncogene, 36(3), 307-318. https://doi.org/10.1038/onc.2016.292.
10. Lopez-Otin, C., & Hunter, T. (2010). The regulatory crosstalk between kinases and proteases in cancer. Nature Reviews Cancer, 10(4), 278-292. https://doi.org/10.1038/nrc2823.
11. Kessenbrock, K., et al. (2010). Matrix metalloproteinases: Regulators of the tumor microenvironment. Cell, 141(1), 52-67. https://doi.org/10.1016/j.cell.2010.03.015.











