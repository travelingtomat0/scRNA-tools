# This program conducts the cell assignments for each subdirectory

# To keep the command running, according to
# https://medium.com/@arnab.k/how-to-keep-processes-running-after-ending-ssh-session-c836010b26a3
# use ssh, then screen, then can disconnect with ctrl-A & ctrl-D and later can resume session with screen -r

# TODO:	Generate some kind of summary file (for each subdirectory): how many cells were analysed?
#	How many cells were removed? How many ensembles? etc. pp.

# TODO: Add the cell-cycle regression?

# TODO: Changed start to 4, since MAHACEB had empty file. Change this back when finishing file!

# TODO: In 38th directory (MYGIFUD) the counts file is corrupt! Wrong format or some other problem.

library(umap)
library(hdf5r)
library(rhdf5)
library(ggplot2)
library(monocle)
library(garnett)
library(org.Hs.eg.db)
library(Seurat)

# Build some base global variables

# The assumption is that this file is run from the "~/infercnv_annotations" folder
# Which contains all relevant directories and files in the resp. sub-directories.
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
marker_genes <- read.table("/home/tomap/Rambow_Markers_Smaller.txt", sep=' ')
m_names <- read.table("/home/tomap/Rambow_Markers_rows_Smaller.txt", sep = ' ')
marker_file_path <- system.file("extdata", "/home/tomap/new_infercnv_annotations/garnett_markup.txt", package = "garnett")
# Iteration over all tumors in the directory
for(i in 20:length(subdirs)) {
# Garbage collection after last iteration
gc()
print(paste("Begin work in directory", subdirs[i], sep = " "))

exc <- try(dir.create(paste(subdirs[i],"/results", sep = "")))

# read in counts data and previous annotations.
infile <- H5Fopen(paste(subdirs[i],"/count_data.h5", sep = ""))
counts <- data.matrix(infile$raw_counts)
rownames(counts) <- make.unique(infile$gene_attrs$gene_names)
colnames(counts) <- infile$cell_attrs[[1]]
# Filter using nexus-annotations
if(head(file.exists(paste(subdirs[i],"/nexus_annotations.txt", sep = "")),1)) {
	t <- try(read.table(paste(subdirs[i],"/nexus_annotations.txt", sep = ""), sep='\t'))
	if(!inherits(t, "try-error")) {
		counts <- counts[,colnames(counts) %in% t[,1]]
	}
}
# Filter non melanoma using infercnv output
# ATTENTION: This will probably remove a lot of cells!
if(head(file.exists(paste(subdirs[i],"/new_annotations.txt", sep = "")),1)) {
	annot <- try(read.table(paste(subdirs[i],"/new_annotations.txt", sep = ""), sep = ','))
	if(!inherits(annot, "try-error")) {
		counts <- counts[,colnames(counts) %in% annot[,1]]
	}
}

m <- marker_genes
colnames(m) <- c("Pigmented","Invasive","NCSC","SMC")
rownames(m) <- m_names[,1]

print("Finsished reading in files")

#counts <- counts[,order(colnames(counts))]
#pd <- data.frame(integer(ncol(counts)))
#rownames(pd) <- colnames(counts)
#colnames(pd) <- c("NonCancerous")
#fd <- data.frame(integer(nrow(counts)))
#rownames(fd) <- rownames(counts)
#pd <- new("AnnotatedDataFrame", data = pd)
#fd <- new("AnnotatedDataFrame", data = fd)

#cds <- newCellDataSet(as(counts, "dgCMatrix"), phenoData = pd, featureData = fd)
#ncds <- estimateSizeFactors(cds)
c <- new_cell_data_set(counts)
#c <- estimateSizeFactors(c)
marker_check <- check_markers(c, "~/new_infercnv_annotations/garnett_markup.txt", db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")

# This saves the marker check in the results section.
# To save a file that checks marker gene expression for the specific clusters, use the commented code below.
# check_markers(new_cell_data_set(data.matrix(counts[,invasive])), "~/new_infercnv_annotations/garnett_markup.txt", db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
# TODO: Can save a table containing ambiguity, marker score etc. Its all in the marker_check object.
# I.e. use code: write.table(data = marker_check, file = paste(subdirs[i],"/results","/marker_information.txt", sep = ""), sep = ',')
pdf(paste(subdirs[i],"/results","/marker_file.pdf", sep = ""))
plot_markers(marker_check)
dev.off()

pbmc_classifier <- train_cell_classifier(cds = c, marker_file = "~/new_infercnv_annotations/garnett_markup.txt",db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", num_unknown = 50, marker_file_gene_id_type = "SYMBOL")

c <- classify_cells(c, pbmc_classifier, db = org.Hs.eg.db, cluster_extend = FALSE, cds_gene_id_type = "SYMBOL", verbose = TRUE)

if(nrow(m) < ncol(counts)) {
umap_model <- umap(t(x), random_state = 42, transform_seed = 69)
df <- data.frame(umap_model$layout)

print("Begin with data fitting.")
set.seed(128)

# Add unassigned state and conduct same experiments
# print("Add unassigned state.")
fit <- cellassign(exprs_obj = t(x), marker_gene_info = data.matrix(m), s = s, learning_rate = 1e-3, shrinkage = TRUE, verbose = FALSE, max_iter_adam = 1500, max_iter_em = 10, min_delta = 1)
type_ass <- data.frame(fit$cell_type)
rownames(type_ass) <- colnames(x)
colnames(type_ass) <- c("Cell Type")

plot <- ggplot(data = df, mapping = aes(x=X1,y=X2, colour = fit$cell_type)) + geom_point() + scale_colour_discrete(drop=TRUE, limits = levels(c("Pigmented","NCSC","Invasive","SMC", "Unassigned"))) + labs(title = "Transcriptional Program Assignments", caption = "Learning rate given to Cellassign: 1e-3. No cycle removal.")
ggsave(paste(subdirs[i],"/results","/assignment3_plot.png", sep = ""), plot = plot)
write.table(type_ass, file = paste(subdirs[i],"/results","/assignment3_cell_types.txt", sep = ""))
write.table(fit$mle_params$gamma, file = paste(subdirs[i],"/results","/assignment3_probs.txt", sep = ""))
print(paste("Finished work in directory", subdirs[i], sep = " "))

## Cell Cycle Effects removal
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

obj <- CreateSeuratObject(counts = counts)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst")
obj <- ScaleData(obj, features = rownames(obj))

# PC8 and PC10 are split on cell-cycle genes 
# [according to https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html]
obj <- RunPCA(obj, features = VariableFeatures(obj), ndims.print = 6:10, nfeatures.print = 10)

obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
# head(obj[[]])
obj <- RunPCA(obj, features = c(s.genes, g2m.genes))


obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj), nfeatures.print = 10)
obj <- RunPCA(obj, features = c(s.genes, g2m.genes))

# Visualize cell cycle effects using counts data
# umap_model <- umap(t(counts), random_state = 42, transform_seed = 69)
# df <- data.frame(umap_model$layout)
# ggplot(data = df, mapping = aes(x=X1,y=X2, colour = obj$Phase)) + geom_point()
# ----------------------

# To use regressed cell cycle data:
# !!! Counts and Object should have same number of genes / cells (?)
counts <- data.matrix(GetAssayData(object = obj, assay.type = "RNA", slot = "counts"))
x <- subset(counts, rownames(counts) %in% m_names[,1])
x <- x[order(rownames(x)),]
x <- x[which(rowSums(x) != 0),]
x <- x[,which(colSums(x) != 0)]

s <- calculateSumFactors(counts)
s <- data.frame(s)
rownames(s) <- colnames(counts)
s <- s[rownames(s) %in% colnames(x),]

m <- m[rownames(m) %in% rownames(x),]

fit <- cellassign(exprs_obj = t(x), marker_gene_info = data.matrix(m), s = s, learning_rate = 1e-3, shrinkage = TRUE, verbose = FALSE, max_iter_adam = 1500, max_iter_em = 10, min_delta = 1)
type_ass <- data.frame(fit$cell_type)
rownames(type_ass) <- colnames(x)
colnames(type_ass) <- c("Cell Type")

plot <- ggplot(data = df, mapping = aes(x=X1,y=X2, colour = fit$cell_type)) + geom_point() + scale_colour_discrete(drop=TRUE, limits = levels(c("Pigmented","Invasive","NCSC","SMC", "Unassigned"))) + labs(title = "Transcriptional Program Assignments", caption = "Learning rate given to Cellassign: 1e-3. With cell cycle removal.")
ggsave(paste(subdirs[i],"/results","/assignment4_plot.png", sep = ""), plot = plot)
write.table(type_ass, file = paste(subdirs[i],"/results","/assignment4_cell_types.txt", sep = ""))
write.table(fit$mle_params$gamma, file = paste(subdirs[i],"/results","/assignment4_probs.txt", sep = ""))
print(paste("Finished work in directory", subdirs[i], sep = " "))

} else {
	write.table(c("Sample from directory", subdirs[i], "contains too little number of samples."), file = paste(subdirs[i], "/warning.txt",sep = ""), sep = " ")
}
}
