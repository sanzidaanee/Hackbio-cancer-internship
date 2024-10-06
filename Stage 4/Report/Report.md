
# Differential Gene Expression Profiling and Unsupervised Clustering to Predict Gliomas Prognosis 



Authors (@slack): Sanzida Akhter Anee (@Sanzida), Sk Arif (@arif_shaikh), Nada Ghozlan (@Nad1), Mennatallah Mohamed Ebrahim Mahmoud (@Mennatallah), Stéphie Raveloson (@StephieRav), Chidimma Nwaku (@Mma),  Zaka Ullah (@Zaka),  Idahosa Clinton (@doc_idahosa)



## Link

#### Github repo 
(https://github.com/sanzidaanee/Hackbio-cancer-internship/tree/main/Stage%204)

### Github code link

DE_Analysis R Code: (https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/Stage%204/Code/DEA.Rmd)

ML Python Code: 

#### Link of dataset

(https://github.com/sanzidaanee/Hackbio-cancer-internship/tree/main/Stage%204/Data)

## Introduction

Gliomas are brain tumors that originate from glial cells, which represent 80% of all primary malignant brain tumors in the central nervous system (CNS) [1].

There are different grades of gliomas, ranging from low-grade (Grade I-II), which are typically slower-growing, than high-grade gliomas (Grade III-IV), which are more aggressive and harder to treat [2]. Mutations occurring in the isocitrate dehydrogenase (IDH) genes (IDH1 and IDH2) identify a subset of glioblastomas with a favorable prognosis [3]. 


## Methods
### Gene expression profiles data
RNA- seq data of low-grade gliomas (LGG)  along with clinical information IDH mutations statuses IDH1 and IDH2 with 516 tumor samples from cancer patients were downloaded from The Cancer Genome Atlas (TCGA) [4]. The selected datasets in total 60660 genes across 534 patient samples of which 460 had mutant IDH and 34 had wild type IDH. 

### Data processing
Through the quantile normalization method by using gene length in edgeR package [5] and log2 transformation of expression matrixes, all datasets were normalized and filtered.

### Differential expression analysis 
Using edgeR package [5], differentially expressed genes were discovered. We used p-value = 0.01 and FDR < 0.01 as  a cutoff value and log  2 fold change (logFC) cutoff value is 1, where >1 for upregulated  and < -1 for downregulated genes.

### Enrichment of DEGs 
By using biomaRt package [6] for gene annotation and TCGAbiolinks package [7] for pathway enrichment analysis at gene ontology biological process level.

### Machine Learning Model

## Results 

### Differential expression analysis

 - To decipher the heterogeneity within IDH-mutant for low-grade gliomas, we filtered 34539 genes with 419 wild types and 94 mutant samples from 60660 genes.

 - Then from 34539 filtered genes, after DE analysis, we reported 5916 genes.

  - We screened all the DEGs using  using log2FC = 0.01 and p-value = 1 as the threshold and showed them as volcano plots (Fig. 1).


    ![Volcano_plot](https://github.com/user-attachments/assets/1573fed5-a214-4d91-a462-41e785319af3)



 Figure 1: Volcano plot of low-grade gliomas (LGG) dataset visualizes gene expression data. X-axis shows log2-fold change and Y-axis shows adjusted p-value. Red points indicate upregulated and blue dots indicate downregulated genes.

 - Following overlapping, we identified 1681 upregulated and 4235 downregulated common genes
 - The results of DE with 5916 genes were  plotted as a heatmap (Fig).

 - Results from heatmap show the clustering of upregulated genes is expressed at a higher level in mutant IDH than wildtype and this result is consistent with the other researcher’s findings (Fig. 2).

   ![heatmap_output](https://github.com/user-attachments/assets/782ec6fe-4585-4e67-99a7-0cc18bb63851)


 Figure 2: Hierarchical clustering heat map of 5916 DEGs from IDH mutant and wild type for low-grade gliomas (LGG) samples, Red and green indicate high and low expression genes, respectively.		

 ### Enrichment analysis

 - Upregulated genes  were enriched in synaptic transmission, cell-cell signaling, transmission of nerve impulse, homophilic cell adhesion and cell-cell adhesion might  tend to have a better prognosis compared to more aggressive gliomas (Fig. 3)

- Downregulated genes were enriched in anterior or posterior pattern formation, skeletal system development, regionalization, pattern specification process and embryonic morphogenesis indicate that these genes have better prognosis than aggressive gliomas (Fig. 4).


![Upregulated_EA](https://github.com/user-attachments/assets/18cdcad1-9e9d-4ceb-aa31-3ac948c0f2a4)






Figure 3: Functional enrichment analysis of upregulated genes for low-grade gliomas (LGG).

![Downregulated_genes_EA](https://github.com/user-attachments/assets/c4a9bac8-97f1-4857-bde0-63c50fdd3851)




Figure 4: Functional enrichment analysis of downregulated genes for low-grade gliomas (LGG).


## Machine learning models



## Conclusion 

  - The results confirm IDH status as the major determinant of the molecular footprints of low grade gliomas. 
 - Mutant LGGs, particularly those with IDH mutations, shows unique gene expression profiles compared to wildtype tumors
 - Enriched pathways suggest that IDH-mutant gliomas might associated with more aggressive, invasive behavior of tumor than wild type




## References

1. Ahir, B. K., Engelhard, H. H., & Lakka, S. S. (2020). Tumor development and angiogenesis in adult brain tumors: glioblastoma. Molecular neurobiology, 57, 2461-2478.

2. What Makes a Brain Tumor High-Grade or Low-Grade? David Reardon, MD. April 3, 2018. https://blog.dana-farber.org/insight/2018/04/makes-brain-tumor-high-grade-low-grade/.
3. Nakhate, V., Lasica, A. B., & Wen, P. Y. (2024). The Role of Mutant IDH Inhibitors in the Treatment of Glioma. Current Neurology and Neuroscience Reports, 1-13.
4. The Cancer Genome Atlas. (2015). Data Portal. National Cancer Institute. https://tcga-data.nci.nih.gov/tcga/tcgaHome2.jsp
5. McCarthy DJ, Chen Y, Smyth GK (2012). "Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation." Nucleic Acids Research, 40(10), 4288-4297. DOI:10.1093/nar/gks042.
6. Durinck, S., Spellman, P.T., Birney, E., & Huber, W. (2009). Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Bioinformatics, 25(18), 2410-2411. https://doi.org/10.1093/bioinformatics/btp515.
7. Colaprico, A., Silva, T. C., Olsen, C., Garofano, L., Cava, C., Garolini, D., ... & Ceccarelli, M. (2016). TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Research, 44(8), e71. https://doi.org/10.1093/nar/gkv1507.













