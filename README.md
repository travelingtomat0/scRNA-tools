# [scRNA-tools]

This repository contains self-programmed tools for scRNA-seq experiments. Also, below is a list of tools for scRNA experiments listed by their use-cases.
For validation purposes there are some useful metrics that I have also listed and described below.

## Contents and Summary of Toolbox

## Cell-Type Annotation
Cell type annotation can be divided into supervised, semi-supervised (prior-knowledge) and unsupervised methods.

Supervised:
- CNN's / ANN's???

Semi-Supervised:
- Cellassign --> performed consistently well on TuPro data set.
- SCINA --> Important note: while cellassign takes a binary marker matrix as input, SCINA has its own preprocessing function which takes a csv marker-file and converts it to a marker gene list.
- Garnett --> not used, due to issues with classification. Also contains a tool to evaluate marker genes.
- scSorter --> unpublished method. GitHub https://github.com/hyguo2/scSorter contains the code, https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html contains a tutorial.

Unsupervised:
- dimred with PCR, UMAP
- RaceID, SC3, CIDR
- RCA
- SIMLR --> clustering method.

## Validation Metrics
Given some gold-standard, we can always use standard validation techniques and metrics. In absence of labeled data, cluster-validation techniques may be useful for validation purposes of cluster assignments (i.e. cell type annotations). Below I have listed the validation metrics that were deployed as part of my Bachelor Thesis.
*For measuring similarity between cells based on scRNA-seq, different papers suggest that using correlation as a measure will give the best results.* 

- *Intracluster compactness* measures the degree of intraclass similarity. The higher the compactness, the greater the intraclass similarity of each cell type. 
- *Intercluster Complexity:* measures the degree of interclass similarity. The lower the complexity, the the smaller is the interclass similarity.
-> High intracluster compactness and low intercluster complexity imply a good assignment performance on the dataset.

- Intersample intracluster compactness can be used to measure the overall assignment performance for our specific task. If the cells assigned to a state are overall similar regardless of the sample they stem from, our assignments are more likely to be correct.

- For actual validation, we need to look at overexpression of type-specific marker genes for assigned clusters.

## Paper-Collection
- Impact of similarity metrics on single-cell RNA-seq data clustering
- Autoencoder-based cluster ensembles for single-cell RNA-seq data analysis
- Learning for single-cell assignment
- (Bayesian Correlation is a robust similarity measure for single-cell RNA-seq data)
- Integrating Deep Supervised, Self-Supervised andUnsupervised Learning for Single-Cell RNA-seqClustering and Annotation
- Cell Type Identification from Single-CellTranscriptomic Data via Semi-supervised Learning

## Framework-Collection
- Seurat (for normalization, analysis of feature expression  etc.)
- scran (for size factor estimation)
