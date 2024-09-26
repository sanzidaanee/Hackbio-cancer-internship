
# Identify Potential Biomarkers for Cancer Detection by Analyzing Differential Expression and Machine Learning Models


Authors (@slack): Sanzida Akhter Anee (@Sanzida), Sk Arif (@arif_shaikh), Nada Ghozlan (@Nad1), Mennatallah Mohamed Ebrahim Mahmoud (@Mennatallah), Stéphie Raveloson (@StephieRav), Chidimma Nwaku (@Mma),  Zaka Ullah (@Zaka),  Idahosa Clinton (@doc_idahosa)



## Link

#### Github repo : (https://github.com/sanzidaanee/Hackbio-cancer-internship/tree/main/Stage%203)
#### Github code link
DE_Analysis R Code: (https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/Stage%203/Code/DEG.Rmd)

ML Python Code: (https://github.com/sanzidaanee/Hackbio-cancer-internship/blob/main/Stage%203/Code/ML%20lung%20Cancer%20.ipynb)

#### Link of dataset
(https://github.com/sanzidaanee/Hackbio-cancer-internship/tree/main/Stage%203/Data)

## Introduction

Lung adenocarcinoma represents a histopathological subtype of non-small cell lung cancer (NSCLC) and has various forms, such as solid, acinar, lepidic, papillary, and micropapillary subtypes [1]. This subtype comprises nearly half of all lung cancer cases [2]. It originates from the glandular epithelial cells that line the alveoli of the lungs and is characterized by the formation of glandular structures or mucin production[3]. 

## Dataset

The lung adenocarcinoma (LUAD) is one of the 3 lung tumors studied by The Cancer Genomic Atlas (TCGA) through its lung cancer research. By June 2015, 521 LUAD samples had been analyzed, and the data had been uploaded to the TCGA data portal. A subset of 230 tumors was the subject of the TCGA lung cancer report for LUAD [4]. 


### Download data

TCGAbiolinks is an R/Bioconductor package designed to facilitate the retrieval, analysis, and integration of data from The Cancer Genome Atlas (TCGA) and other Genomic Data Commons (GDC) resources.

## Data Preprocessing

### Step 1: Create metadata

Metadata is important for exploring and making meaningful insights from the large dataset. Choose metadata based on race, tumor type, sample type and id.


### Step 2: Downsize data set

20 select for primary tumor and 20 for Solid Tissue Normal data

### Step 3: Data normalization and filtering


#### Data Normalization

- Normalization in RNA-Seq data is used to adjust raw read counts to account for various biases and ensure that gene expression levels can be accurately compared across different samples.


#### Filtering

- Filtering is important to filter out low-expressed genes from a normalized gene expression matrix. The quantile filtering method works by calculating a specified quantile threshold for the gene expression data and removing genes with expression values below this threshold.

- In this case, quantile cut off value 0.25 means that the bottom 25% of genes with the lowest expression values across all samples will be removed as they are likely to be noise or irrelevant to the analysis.

## Methods for Biomarker Discovery


### Differential Expression Analysis


Differential expression analysis is a fundamental technique used to identify changes in gene expression between different biological conditions like healthy or diseased situations. By analyzing the RNA-seq data can identify the up and down regulated genes that respond to specific conditions or treatments can be used as a potential biomarker for early detection of cancer [5]. 


### Functional Enrichment Analysis

Enrichment analysis is used to identify biological pathways, gene sets or functional categories in a given set of genes that are associated with specific types of cancer visualizing by barplot focusing on biological processes, cellular components, molecular functions, and pathways associated with genes. 


## Results and Visualization

### Differentially Expressed Genes  (DEGs)


- DEGs were identified using a t-statistical test, setting a significance level of 1%. The log2 fold change (FC) was calculated with cutoff 1 and >1 for upregulation and <   -1 for downregulation. 

- From this analysis, 3277 upregulated genes and 6357 downregulated genes were extracted from a raw extracted dataset. 
- A volcano plot was created to visualize these results, showing how highly differentially  genes were expressed based on fold changes and p values. So, the volcano plot would show genes up-regulated in lung adenocarcinoma on the right side (e.g., cancer-associated genes), while those down-regulated in the cancer state appear on the left side.
  
![Volcano Plot](https://github.com/user-attachments/assets/17e6767b-6904-4545-a278-828a3e3b72f4)







Figure  1: Volcano plot of lung adenocarcinoma dataset visualize gene expression data. X-axis shows log2-fold change and Y-axis shows adjusted p-value. Red points indicate upregulated and blue dots indicate downregulated genes








- The results of  heatmap reveal the tumor and normal solid tissue samples clustering separately on both sides of the heatmap suggest a strong differential gene expression between two conditions (tumor and normal solid tissue). Genes that are expressed at higher levels in tumor tissues appear red in the tumor column and blue in normal tissue and these genes  might be involved in cancer progression or oncogenic processes. 

- On the other hand,  genes that are expressed at lower levels in tumor tissues will appear blue in tumor columns and red in normal tissue columns may be a tumor suppressor gene.


So, identifying the key genes that are consistently differentially expressed between tumor and normal tissues, can be used as potential biomarkers for diagnosis or prognosis.


![Heatmap](https://github.com/user-attachments/assets/1a25877d-dec4-4e2c-a9c7-15c5fc11906b)





Figure 2: Heatmap showing clusters on both genes and samples



## Functional Enrichment Analysis

- Functional enrichment analysis was conducted on differentially expressed genes from lung cancer samples, revealing significant enrichment of several key biological processes. 

- In upregulated genes, the top three GO biological process are included behavior,  cell adhesion and biological adhesion and their upregulation is highly relevant in lung cancer because these mechanisms enable tumor cells to become invasive, migrate to distant sites, and ultimately lead to metastasis which is a major cause of lung cancer mortality


- On the other hand, in cases of down regulated genes, the  top biological processes are cell-cell signaling, nuclear division, and mitosis in lung adenocarcinoma (LUAD) indicates disruptions in key regulatory pathways that control cellular communication, cell division, and proliferation. 


![Upregulated_genes_EA](https://github.com/user-attachments/assets/fef42a14-417a-45a3-b8fc-ee42c32bd3c2)




Figure 3: Functional enrichment analysis of upregulated genes for lung adenocarcinoma


![Downregulated_genes_EA](https://github.com/user-attachments/assets/fb3d3ccb-c456-41da-bddf-ac7dae77f3c3)



Figure 4: Functional enrichment analysis of downregulated genes for lung adenocarcinoma


## Machine Learning Models

Machine learning can be used in cancer research to distinguish between tumor and solid normal tissue. ML models can help in early detection of cancer by comparing tumor and normal tissues, and identify novel biomarkers for cancer diagnosis. 


The Random Forest model is a widely used machine learning algorithm for classification of gene expression data. This method constructs multiple decision trees during training and outputs the mode of the classes (classification) or mean prediction (regression) that improve accuracy. Random feature selection is considered for splitting the dataset.

 - Split the dataset into training and testing sets (e.g., 80% for training, 20% for testing).
 - Check for duplicate columns both in the training and test data
 - Transposed normalized dataset 
 - Renamed the first column to “barcode” and drop the first row
 - We used Pandas and RandomForest Classifier libraries respectively to manipulate and prepare the dataset and train the model using Random 
    Forest.
 - Train the Random Forest Model: randomly select subset of the data and train  to classify whether a tissue sample is primary or solid 
   based on the input features.
 - Test the Model: After training, the model is tested to evaluate its performance.



## Conclusion 

- DEG analysis revealed a set of upregulated and downregulated genes in lung adenocarcinoma (LUAD), providing insight into the molecular mechanisms underlying cancer progression

- The subsequent enrichment analysis of these genes offers valuable information on the biological processes, molecular functions, and pathways that are potentially disrupted during tumor development

- ML models can be trained on gene expression data to  identify key biomarkers for lung cancer that are linked to prognosis and  targeted therapy



## References

1. Liang, Y., Xie, Y., Yu, H., Zhu, W., Yin, C., Zhang, X., & Dong, Z. (2023). Clinical significance of TMEM229A Q200del mutation in lung adenocarcinoma.
2. Shiba‐Ishii, A. (2021). Significance of stratifin in early progression of lung adenocarcinoma and its potential therapeutic relevance. Pathology international, 71(10), 655-665.
3. Herbst, R. S., Morgensztern, D., and Boshoff, C. "The biology and management of non-small cell lung cancer." Nature, vol. 553, 2018, pp. 446-454.
4. The Cancer Genome Atlas. (2015). Data Portal. National Cancer Institute. https://tcga-data.nci.nih.gov/tcga/tcgaHome2.jsp

5. Wang, Z., Gerstein, M., and Snyder, M. "RNA-Seq: a revolutionary tool for transcriptomics." Nature Reviews Genetics, vol. 10, no. 1, 2009, pp. 57-63.








