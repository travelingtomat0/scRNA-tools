# Summary of the Code Snippets Directory

- [get_abs_cnv.py] takes a cnv expression file (either one where levels are 1 - 6 (1 resembling high loss of gene copy numbers, 6 high copynumber overexpression of gene) or one of distribution centered around one). The value around which is centered needs to be taken into account in th code ;-)
  The code returns the mean copynumber variation either for all genes or only for marker genes. This also needs to be edited in the code depending on what is to be acheived.

- [full_analysis.R] generates cellassign outputs and visualizes a umap plot, using Nexus and infercnv to remove non melanocytic / mesenchymial cells. Can be modified to run with different states by editing the marker matrix.

- [full_validation.R] generates validations based on validation metrics (see README.md).

- [visualize_all_tumors.R] visualizes all tumor samples in one umap plot. It is possible to either run umap on all features, or run umap only on the marker genes measured for all samples (which can be less than the actual number of marker genes originally specified).

- [ensemble_autoenc_code.R] generates an autoencoder generated encoding, on which k-means (in an ensemble clustering setting) is then run. It's possible to specify other clustering methods!

- [subclonal_cnv_analysis.R] takes an infercnv dendrogram as an input (and raw counts, which can be omitted if script is changed such that previously calculated umap-plot is used). It then removes non melanocytic cells according to Nexus and Infercnv annotations, extracts the subclonal structures (up to level k, which can be specified in the call of subtree_split(tree, 1, k)) and then visualizes the plot based on the subclonal structures found by infercnv.
