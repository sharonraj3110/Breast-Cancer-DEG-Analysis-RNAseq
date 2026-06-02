&#x20;**Breast Cancer RNA-seq Differential Expression Analysis**



* Project Overview



This project analyzes RNA-seq gene expression data from breast cancer samples (GEO dataset: GSE183947) to identify differentially expressed genes between tumor and normal tissues using DESeq2 in R.



The goal is to detect potential biomarkers and understand gene expression changes associated with cancer progression.



* &#x20;Objectives

Perform differential gene expression analysis

Identify significantly upregulated and downregulated genes

Visualize expression patterns using plots

Interpret biological relevance of results



* &#x20;Dataset

Source: NCBI Gene Expression Omnibus (GEO)

Accession: GSE183947

Type: RNA-seq

Groups: Tumor vs Normal breast tissue



* &#x20;Tools \& Technologies

R Programming

DESeq2 (Differential Expression Analysis)

GEOquery (Data acquisition)

ggplot2 (Visualization)

pheatmap (Heatmap generation)



* &#x20;Methodology

Data acquisition from GEO database

Preprocessing and normalization of raw counts

Differential expression analysis using DESeq2

Identification of significant genes (padj < 0.05, |log2FC| > 1)

Data visualization using:

Volcano plot

Heatmap

PCA plot



* Results

Identified significant upregulated genes in tumor samples

Identified significant downregulated genes in tumor samples

Clear separation between tumor and normal samples in PCA analysis

Expression patterns visualized using heatmap clustering



* &#x20;Repository Structure

scripts/ → R analysis pipeline

data/ → gene expression results

figures/ → plots and visualizations

README.md → project documentation



* &#x20;Biological Insight

The analysis highlights dysregulation of genes involved in key cancer-related processes such as cell cycle regulation, DNA repair, and apoptosis, indicating their potential role in breast cancer progression.



* Key Skills Demonstrated

RNA-seq data analysis

Statistical modeling in R

Bioinformatics workflow design

Data visualization

Biological interpretation of gene expression data

