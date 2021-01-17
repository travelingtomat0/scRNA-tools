# CNV tumor subclonal structure analysis
library(hdf5r)
library(rhdf5)
library(ggplot2)
library(umap)
library(ape)

# Reading: http://www.phytools.org/eqg/Exercise_3.2/

# Take the dendrogram, which is given in Newick Format. Process the tree,
# such that tips (leafs) are removed from the tree.
global_count <- 1

tree_split <- function(tree, curr_level, end_level, k) {
	if(nrow(tree$edge) == 1) {
		return(data.frame(c()))
	}
	if(curr_level == end_level) {
		assignments <- data.frame(integer(length(tree$tip.label)) + global_count)
		rownames(assignments) <- tree$tip.label
		tmp <- global_count + 1
		assign("global_count", tmp, envir = .GlobalEnv)
		colnames(assignments) <- c("Family")
		return(assignments)
	}
	# find the subtree roots of the subtrees
	root_node <- tree$edge[1,1]
	child_roots <- c(tree$edge[1,2])
	for(i in 2:nrow(tree$edge)) {
		if(tree$edge[i,1] == root_node) {
			child_roots <- c(child_roots, tree$edge[i,2])
		}
	}
	# return list of all subtrees that have a child of the current root as their root.
	res <- data.frame(c())
	for(r in child_roots) {
		if(r != 1) {
		res <- rbind(res, tree_split(extract.clade(tree, r), curr_level + 1, end_level))
		}
	}
	return(res)
}

infile <- H5Fopen("count_data.h5")
counts <- infile$raw_counts
rownames(counts) <- make.unique(infile$gene_attrs$gene_names)
colnames(counts) <- infile$cell_attrs[[1]]

t <- read.tree("infercnv.observations_dendrogram.txt")
n <- read.table("nexus_annotations.txt", sep='\t')
annot <- read.table("new_annotations.txt", sep='\t')

for(i in setdiff(colnames(counts), c(n[,1],annot[,1]))) {
	t <- drop.tip(t, i)
}

counts <- counts[,colnames(counts) %in% n[,1]]
counts <- counts[,colnames(counts) %in% annot[,1]]
counts <- data.matrix(counts)

# Now analyse tree structure to get the  subclonal structure.

subtrees <- tree_split(t, 1, 2)
# TODO: remove barcodes not found in counts:
subtrees <- subset(subtrees, rownames(subtrees) %in% colnames(counts))
# TODO: use subtrees[order(rownames(subtrees)),] when calling ggplot, so that assignments are colored by assignments correctly.

# for each of the subclonal structures, we can get a list of the tips ("leafs"), which are labeled with the cell barcodes.
# we then generate a dataframe, which contains all the cells with their resp. families as number labels.

# Color umap by the family


#umap_plot <- read.table("umap_layout")
m_names <- read.table("/home/tomap/Rambow_Markers_rows_Smaller.txt", sep = ' ')
x <- subset(counts, rownames(counts) %in% m_names[,1])
umap_model <- umap(t(x), random_state = 42, transform_seed = 69)
# Probably a good idea to plot non reduced plot as well!
# umap_model <- umap(t(counts), random_state = 42, transform_seed = 69)
df <- data.frame(umap_model$layout)
subclonal_structure <- factor(subtrees[order(rownames(subtrees)),])
plot <- ggplot(data = df, mapping = aes(x=X1,y=X2, colour = subclonal_structure)) + geom_point() + scale_colour_discrete(drop=FALSE, limits = levels(unique(subtrees[,"Family"]))) + labs(title = "Subclonal Structure", caption = "Dendrogram taken from infercnv output. Structures found at level 3 of the hierarchy.")
