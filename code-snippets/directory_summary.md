# Summary of the Code Snippets Directory

- [get_abs_cnv.py] takes a cnv expression file (either one where levels are 1 - 6 (1 resembling high loss of gene copy numbers, 6 high copynumber overexpression of gene) or one of distribution centered around one). The value around which is centered needs to be taken into account in th code ;-)
  The code returns the mean copynumber variation either for all genes or only for marker genes. This also needs to be edited in the code depending on what is to be acheived.

- [full_analysis.R]

- [full_validation.R]

- [visualize_all_tumors.R]

- [ensemble_autoenc_code.R]

- [subclonal_cnv_analysis.R] takes an infercnv dendrogram as an input (and raw counts, which can be omitted if script is changed such that previously calculated umap-plot is used). It then removes non melanocytic cells according to Nexus and Infercnv annotations, extracts the subclonal structures (up to level k, which can be specified in the call of subtree_split(tree, 1, k)) and then visualizes the plot based on the subclonal structures found by infercnv.
