# [scRNA-tools]

This repository contains self-programmed tools for scRNA-seq experiments. Also, below is a list of tools for scRNA experiments listed by their use-cases.
For validation purposes there are some useful metrics that I have also listed and described below.

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
Given some gold-standard, we can always use standard cluster-validation techniques for validation purposes of cluster assignments (i.e. cell type annotations). For other cases (absence of labels or equivalent data) I have developed other validation metrics (as part of my Bachelor Thesis).

- 

# Paper-Collection
- Impact of similarity metrics on single-cell RNA-seq data clustering
- Autoencoder-based cluster ensembles for single-cell RNA-seq data analysis
- Learning for single-cell assignment
- (Bayesian Correlation is a robust similarity measure for single-cell RNA-seq data)
