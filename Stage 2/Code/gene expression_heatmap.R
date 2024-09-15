# Part1: Data pre-processing


## Install and load Required Packages


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('biomaRt')
install.packages("gplots")
install.packages("dplyr")
install.packages("ggplot2")




## Load required libraries



llibrary(biomaRt)
library(gplots)
library(dplyr)
library(ggplot2)




## Download data


### Gene Expression Dataset

- Glioblastoma gene expression dataset with 500+ differentially expressed genes under different conditions.

- Glioblastoma, also known as glioblastoma multiforme (GBM), is an aggressive type of cancer that occurs in the brain or spinal cord. The data represent the samples of cancer tissue from the patients with gene expression data for many different types of genes.



data <- read.csv ("/Users/sanzidaakhteranee/Documents/Internship_2024/HackBio/sTAGE 2/glioblastoma.csv", header = TRUE, row.names = 1 )


## Data structure


head(data)
View(data)
colnames(data)       
nrow(data)             
ncol(data)




# Part 2: Heatmap Generation

-Creating two color variants of the same heatmap by using diverging and sequential color palettes in gene expression data, with heatmap.2() function from the gplots package in R. The choice of color palette is important for interpreting gene expression data effectively to distinguish between different patterns in gene regulation.

Two color palettes:
  
  -Diverging color palette (e.g., green to brown): represent both upregulation and downregulation in gene expression
-Sequential color palette (e.g., white to blue) : highlight the intensity of gene expression, irrespective of the direction of regulation

## Generate heatmap


heatmap.2(as.matrix(data), trace = 'none')





## Scaling


heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE)


## Heatmap with diverging color palettes


heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE,
          col=hcl.colors(100, palette = 'green-brown'))



## Heatmap with sequential color palettes


heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE,
          col=hcl.colors(100, palette = 'Blues3'))



# Part3: Heatmap Clustering

The heatmap can be combined with clustering methods which represent group genes with or without samples together based on the similarity of their gene expression pattern. This is useful for identifying genes that are commonly regulated, or biological signatures associated with a particular disease condition

## Heatmap with Clustering of Genes (Rows Only)


heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'row', 
          Colv = FALSE, Rowv = TRUE,
          col=hcl.colors(100, palette = 'green-brown'))



## Heatmap with Clustering of Samples (Columns Only)


heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'col', 
          Colv = TRUE, Rowv = FALSE,
          col=hcl.colors(100, palette = 'green-brown'))




## Heatmap with Clustering both Genes and Sample Together


heatmap.2(as.matrix(data), trace = 'none', 
          scale='row', dendrogram = 'both', 
          Colv = TRUE, Rowv = TRUE,
          col=hcl.colors(100, palette = 'green-brown'))



# Part 4: Fold Change and P Value Calculation


## Selecting groups by index positions


group1<- c(1, 2,3,4,5)
group2<- c(6, 7,8,9,10)



## Group 1 & 2 from data


group1_data <- data[, group1]
group2_data <- data[, group2]



## Calculate group mean


group1_mean <- rowMeans(group1_data)
group2_mean <- rowMeans(group2_data)



## Calculate fold change


fold_change <- (group2_mean)-(group1_mean)/group1_mean
logFC <- log2(fold_change)


`

## Calculate P value



pvalues <- apply(data, 1, function(row) {
  t.test(row[1:5], row[6:10])$p.value
})





## Add fold change and p-values to the dataset



results <- data.frame(Gene = rownames(data), logFC, pvalues)




# Part 5: Subset genes


## Set cutoff thresholds for fold change and p-value



logFC_up_cutoff <-  1   # Log2 fold change > 1 for upregulation
logFC_down_cutoff <- -1 # Log2 fold change < -1 for downregulation
pvalue_cutoff <- 0.05



## Subset significantly upregulated genes



upregulated_genes <- results %>%
  filter(logFC > logFC_up_cutoff & pvalues < pvalue_cutoff)




## Subset significantly downregulated genes


# Subset significantly downregulated genes
downregulated_genes <- results %>%
  filter(logFC < logFC_down_cutoff & pvalues < pvalue_cutoff)



## Print the results


print(upregulated_genes)
print(downregulated_genes)




## View the first few rows of the upregulated and downregulated genes



head(upregulated_genes)
head(downregulated_genes)





## Save upregulated genes to CSV


write.csv(upregulated_genes, "upregulated_genes.csv")



## Save downregulated genes to CSV


write.csv(downregulated_genes, "downregulated_genes.csv")



## Convert Ensembl IDs to gene id (upregulated genes) using biomaRt

*BioMart is a powerful, web-based data management and analysis tool that integrates data from multiple biological databases, including Ensembl (genome data) and others*
  
  
  ### Connect to the Ensembl BioMart database for human genes
  
  
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)



### List of Ensembl IDs from upregulated genes data


ensembl_ids <- c("ENSG00000241945", "ENSG00000279104", "ENSG00000231107", 
                 "ENSG00000254092", "ENSG00000172236", "ENSG00000197253", 
                 "ENSG00000172116", "ENSG00000162598", "ENSG00000256193", 
                 "ENSG00000160183")



###Retrieve gene symbols for the Ensembl IDs


gene_conversion <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                         filters = "ensembl_gene_id",
                         values = ensembl_ids,
                         mart = ensembl)



### View the results


print(gene_conversion)





#Part 6: Functional Enrichment Analysis

Functional enrichment analysis is a process to identify biological pathways, functions, or processes that are over-represented in a given set of genes that provide information into the biological mechanisms underlying a specific condition, disease progression or response to treatment.


We used the ShinyGO (Version: v0.741) tool with the GO biological process and P-value cutoff (FDR) is 0.5 to visualize the top 10 pathways. From here we choose to describe the top 3 enriched pathways based on biological processes.


## Pathway Visualization in R


### Upload data

-Upregulated data from functional enrichment analysis


enrichment_pathway <- read.csv("/Users/sanzidaakhteranee/Documents/Internship_2024/HackBio/sTAGE 2/pathway_R/enrichment .csv")





# Input data from the table in the image
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

