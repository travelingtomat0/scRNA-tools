# This program conducts the cell assignments for each subdirectory

# To keep the command running, according to
# https://medium.com/@arnab.k/how-to-keep-processes-running-after-ending-ssh-session-c836010b26a3
# use ssh, then screen, then can disconnect with ctrl-A & ctrl-D and later can resume session with screen -r

# TODO:	Generate some kind of summary file (for each subdirectory): how many cells were analysed?
#	How many cells were removed? How many ensembles? etc. pp.

# TODO: Add the cell-cycle regression?

# TODO: Changed start to 4, since MAHACEB had empty file. Change this back when finishing file!

# TODO: In 38th directory (MYGIFUD) the counts file is corrupt! Wrong format or some other problem.

# To get a table with freq. of marker gene expression: t <- table(x["PMEL",]) and then ggplot(data = df, mapping = aes(x = Var1, y = Freq)) + geom_point()

# TODO: The drop command in umap means that not expressed levels are not on scale. This should be false for nicer visuals!!

# IMPORTANT: command to copy file structure find . -name '*_cell_types.txt' | cpio -pdm ~/new_annotations_copy
# find . -name 'assignment*_plot.png' | cpio -pdm ~/new_annotations_copy
# find . -name 'assignment*_probs.txt' | cpio -pdm ~/new_annotations_copy

library(umap)
library(hdf5r)
library(rhdf5)
library(ggplot2)
library(scran)
library(cellassign)

# Build some base global variables

# The assumption is that this file is run from the "~/infercnv_annotations" folder
# Which contains all relevant directories and files in the resp. sub-directories.
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
marker_genes <- read.table("/home/tomap/Rambow_Markers_Smaller.txt", sep=' ')
m_names <- read.table("/home/tomap/Rambow_Markers_rows_Smaller.txt", sep = ' ')

# Iteration over all tumors in the directory
for(i in 1:length(subdirs)) {
# Garbage collection after last iteration
gc()
print(paste("Begin work in directory", subdirs[i], sep = " "))

exc <- try(dir.create(paste(subdirs[i],"/results", sep = "")))

# read in counts data and previous annotations.
infile <- H5Fopen(paste(subdirs[i],"/count_data.h5", sep = ""))
counts <- infile$raw_counts
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
if(head(file.exists(paste(subdirs[i],"/annotations.txt", sep = "")),1)) {
	annot <- try(read.table(paste(subdirs[i],"/annotations.txt", sep = ""), sep = ','))
	if(!inherits(annot, "try-error")) {
		counts <- counts[,colnames(counts) %in% annot[,1]]
	}
}

m <- marker_genes
colnames(m) <- c("Pigmented","Invasive","NCSC","SMC")
rownames(m) <- m_names[,1]

print("Finsished reading in files")
x <- subset(counts, rownames(counts) %in% m_names[,1])
x <- x[order(rownames(x)),]
x <- x[which(rowSums(x) != 0),]
x <- x[,which(colSums(x) != 0)]

s <- calculateSumFactors(counts)
s <- data.frame(s)
rownames(s) <- colnames(counts)
s <- s[rownames(s) %in% colnames(x),]

m <- m[rownames(m) %in% rownames(x),]

if(nrow(m) < ncol(x)) {
umap_model <- umap(t(x), random_state = 42, transform_seed = 69)
df <- data.frame(umap_model$layout)

print("Begin with data fitting.")
set.seed(128)
fit <- cellassign(exprs_obj = t(x), marker_gene_info = data.matrix(m), s = s, learning_rate = 1e-3, shrinkage = TRUE, verbose = FALSE, max_iter_adam = 1500, max_iter_em = 10, min_delta = 1)
type_ass <- data.frame(fit$cell_type)
rownames(type_ass) <- colnames(x)
colnames(type_ass) <- c("Cell Type")

plot <- ggplot(data = df, mapping = aes(x=X1,y=X2, colour = fit$cell_type)) + geom_point() + scale_colour_discrete(drop=FALSE, limits = levels(c("Pigmented","NCSC","Invasive","SMC"))) + labs(title = "Transcriptional Program Assignments", caption = "Learning rate given to Cellassign: 1e-3. No cycle removal.")
ggsave(paste(subdirs[i],"/results","/assignment1_plot.png", sep = ""), plot = plot)
write.table(type_ass, file = paste(subdirs[i],"/results","/assignment1_cell_types.txt", sep = ""))
write.table(fit$mle_params$gamma, file = paste(subdirs[i],"/results","/assignment1_probs.txt", sep = ""))

print("Begin with data fitting 2.")

fit <- cellassign(exprs_obj = t(x), marker_gene_info = data.matrix(m), s = s, learning_rate = 1e-2, shrinkage = TRUE, verbose = FALSE, max_iter_adam = 1500, max_iter_em = 10, min_delta = 1)
type_ass <- data.frame(fit$cell_type)
rownames(type_ass) <- colnames(x)
colnames(type_ass) <- c("Cell Type")

plot <- ggplot(data = df, mapping = aes(x=X1,y=X2, colour = fit$cell_type)) + geom_point() + scale_colour_discrete(drop=FALSE, limits = levels(c("Pigmented","NCSC","Invasive","SMC"))) + labs(title = "Transcriptional Program Assignments", caption = "Learning rate given to Cellassign: 1e-2. No cycle removal.")
ggsave(paste(subdirs[i],"/results","/assignment2_plot.png", sep = ""), plot = plot)
write.table(type_ass, file = paste(subdirs[i],"/results","/assignment2_cell_types.txt", sep = ""))
write.table(fit$mle_params$gamma, file = paste(subdirs[i],"/results","/assignment2_probs.txt", sep = ""))

# Add unassigned state and conduct same experiments
print("Add unassigned state.")
m$Unassigned <- integer(nrow(m))
print("Begin with data fitting 3.")
fit <- cellassign(exprs_obj = t(x), marker_gene_info = data.matrix(m), s = s, learning_rate = 1e-3, shrinkage = TRUE, verbose = FALSE, max_iter_adam = 1500, max_iter_em = 10, min_delta = 1)
type_ass <- data.frame(fit$cell_type)
rownames(type_ass) <- colnames(x)
colnames(type_ass) <- c("Cell Type")

plot <- ggplot(data = df, mapping = aes(x=X1,y=X2, colour = fit$cell_type)) + geom_point() + scale_colour_discrete(drop=FALSE, limits = levels(c("Pigmented","NCSC","Invasive","SMC", "Unassigned"))) + labs(title = "Transcriptional Program Assignments", caption = "Learning rate given to Cellassign: 1e-3. No cycle removal.")
ggsave(paste(subdirs[i],"/results","/assignment3_plot.png", sep = ""), plot = plot)
write.table(type_ass, file = paste(subdirs[i],"/results","/assignment3_cell_types.txt", sep = ""))
write.table(fit$mle_params$gamma, file = paste(subdirs[i],"/results","/assignment3_probs.txt", sep = ""))

print("Begin with data fitting 4.")
fit <- cellassign(exprs_obj = t(x), marker_gene_info = data.matrix(m), s = s, learning_rate = 1e-2, shrinkage = TRUE, verbose = FALSE, max_iter_adam = 1500, max_iter_em = 10, min_delta = 1)
type_ass <- data.frame(fit$cell_type)
rownames(type_ass) <- colnames(x)
colnames(type_ass) <- c("Cell Type")

plot <- ggplot(data = df, mapping = aes(x=X1,y=X2, colour = fit$cell_type)) + geom_point() + scale_colour_discrete(drop=FALSE, limits = levels(c("Pigmented","NCSC","Invasive","SMC", "Unassigned"))) + labs(title = "Transcriptional Program Assignments", caption = "Learning rate given to Cellassign: 1e-2. No cycle removal.")
ggsave(paste(subdirs[i],"/results","/assignment4_plot.png", sep = ""), plot = plot)
write.table(type_ass, file = paste(subdirs[i],"/results","/assignment4_cell_types.txt", sep = ""))
write.table(fit$mle_params$gamma, file = paste(subdirs[i],"/results","/assignment4_probs.txt", sep = ""))
print(paste("Finished work in directory", subdirs[i], sep = " "))
} else {
	write.table(c("Sample from directory", subdirs[i], "contains too little number of samples."), file = paste(subdirs[i], "/warning.txt",sep = ""), sep = " ")
}
}
