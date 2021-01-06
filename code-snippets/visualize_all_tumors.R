# This program visualizes all tumors in one umap plot (features: marker gene raw counts, normalized).
# TODO: We want to visualize both: state and patient! Therefore need to store some kind of tuple-colname, i.e. (Tumor-Name, Barcode).

#TODO??? (substr(dir_name,3,9) or (substr(dir_name,3,10)??? Should be nine, but put 10 , so two '_' after sample..

library(umap)
library(hdf5r)
library(rhdf5)
library(ggplot2)
#library(Seurat)

# Build some base global variables and functions

makeColnames <- function(barcodes, dir_name) {
tmp <- c()
for(code in barcodes) {
	tmp <- c(tmp, paste(substr(dir_name,3,10),code, sep = "_"))
}
return(tmp)
}

# The assumption is that this file is run from the "~/infercnv_annotations" folder
# Which contains all relevant directories and files in the resp. sub-directories.
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
m_names <- read.table("/home/tomap/Rambow_Markers_rows_Smaller.txt", sep = ' ')

res_matrix <- integer(nrow(m_names))
res_matrix <- data.frame(res_matrix)
rownames(res_matrix) <- m_names[,1]
colnames(res_matrix) <- c("Remove")
res_matrix <- t(res_matrix)

# Iteration over all tumors in the directory
for(i in 1:36) {
# Garbage collection after last iteration
gc()
print(paste("Collect data from", subdirs[i], sep = " "))

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

print("Finsished reading in files")

x <- subset(counts, rownames(counts) %in% m_names[,1])
x <- x[order(rownames(x)),]
print(ncol(t(x)))
colnames(x) <- makeColnames(colnames(x), subdirs[i])


res_matrix <- rbind(res_matrix[,colnames(res_matrix) %in% colnames(t(x))], t(x[rownames(x) %in% colnames(res_matrix),]))

print("Finsished appending to results-matrix.")
}

res_matrix <- res_matrix[2:nrow(res_matrix),]

print(ncol(res_matrix))

umap_model <- umap(res_matrix, random_state = 42, transform_seed = 69)
df <- data.frame(umap_model$layout)
write.table(df, "umap_all_tumors.txt", sep = " ")
ggplot(data = df, mapping = aes(x=X1,y=X2, colour = df)) + geom_point()

getSampleVec <- function(rnames) {
	temp <- c()	
	for(i in 1:length(rnames)) {
		temp <- c(temp, substr(rnames[i],1,7))
	}
	return(temp)
}

type <- getSampleVec(rownames(t))
ggplot(data = data.frame(t), mapping = aes(x=X1,y=X2, color = type)) + geom_point()

rownameChanger <- function(rnames, dirname) {
	temp <- c()
	for(i in 1:length(rnames)) {
		temp <- c(temp, paste(dirname,rnames[i], sep = "_"))
	}
	return(temp)
}

ass_matrix <- data.frame(c())

for(i in 1:36) {
	if(head(file.exists(paste(subdirs[i],"/results/assignment3_cell_types.txt", sep = "")),1)) {
	ass <- try(read.table(paste(subdirs[i],"/results/assignment3_cell_types.txt", sep = ""), sep= " "))
	if(inherits(ass, "try-error")) {
		print(paste("Directory",subdirs[i],"doesn't contain assignments", sep = " "))
	} else {
		rownames(ass) <- rownameChanger(rownames(ass), substring(subdirs[i], 3, 10))
		ass_matrix <- rbind(ass_matrix, ass)
	}
}
}

