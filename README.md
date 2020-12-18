# [scRNA-tools]

This repository contains self-programmed tools for scRNA-seq experiments. Also, below is a list of tools for scRNA experiments listed by their use-cases.
For validation purposes there are some useful metrics that I have also listed and described below.

# Contents and Summary of Toolbox

# Cell-Type Annotation
Cell type annotation can be divided into supervised, semi-supervised (prior-knowledge) and unsupervised methods.

Supervised:
- CNN's / ANN's???

Semi-Supervised:
- Cellassign
- SCINA
- Garnett

Unsupervised:
- PCR and its variants
- RaceID, SC3, CIDR
- RCA
- SIMLR

# Validation Metrics
Given some gold-standard, we can always use standard validation techniques and metrics. In absence of labeled data, cluster-validation techniques may be useful for validation purposes of cluster assignments (i.e. cell type annotations). Below I have listed the validation metrics that were deployed as part of my Bachelor Thesis.
*For measuring similarity between cells based on scRNA-seq, different papers suggest that using correlation as a measure will give the best results.* 

- *Intracluster compactness* measures the degree of intraclass similarity. The higher the compactness, the greater the intraclass similarity of each cell type. 
- *Intercluster Complexity:* measures the degree of interclass similarity. The lower the complexity, the the smaller is the interclass similarity.
-> High intracluster compactness and low intercluster complexity imply a good assignment performance on the dataset.

- Intersample intracluster compactness can be used to measure the overall assignment performance for our specific task. If the cells assigned to a state are overall similar regardless of the sample they stem from, our assignments are more likely to be correct.

- For actual validation, we need to look at overexpression of type-specific marker genes for assigned clusters.

# Paper-Collection
- Impact of similarity metrics on single-cell RNA-seq data clustering
- Autoencoder-based cluster ensembles for single-cell RNA-seq data analysis
- Learning for single-cell assignment
- (Bayesian Correlation is a robust similarity measure for single-cell RNA-seq data)
