# use the scCESS autoencoder clustering for our purpose
library(clue)
library(hdf5r)
library(rhdf5)
library(keras)
library(parallel)
library(ggplot2)
library(umap)

source(file = "r_autoencoder.R")

infile <- H5Fopen("count_data.h5")
counts <- infile$raw_counts
rownames(counts) <- make.unique(infile$gene_attrs$gene_names)
colnames(counts) <- infile$cell_attrs[[1]]

t <- read.table("nexus_annotations.txt", sep='\t')
counts <- counts[,colnames(counts) %in% t[,1]]
annot <- read.table("new_annotations.txt", sep='\t')
counts <- counts[,colnames(counts) %in% annot[,1]]
counts <- data.matrix(counts)

df <- encode(t(counts), seed = 128, max_random_projection = ncol(counts), encoded_dim = 32, hidden_dims = c(256), learning_rate = 0.001, batch_size = 32, epochs = 100, verbose = 1, scale = FALSE, genes_as_rows = FALSE)

df <- ensemble_cluster(t(counts), seed = 128, cores = 4, cluster_func = function(x) kmeans(x, centers=4), ensemble_sizes = c(1,5), max_random_projection = ncol(counts), encoded_dim = 16, hidden_dims = c(128), learning_rate = 0.001, batch_size = 32, epochs = 100, verbose = 1, scale = FALSE, genes_as_rows = FALSE)

bf <- ensemble_cluster(t(counts), seed = 128, cores = 4, cluster_func = function(x) kmeans(x, centers=4), ensemble_sizes = c(1,5,10), max_random_projection = ncol(counts), encoded_dim = 32, hidden_dims = c(128,512), learning_rate = 0.001, batch_size = 32, epochs = 100, verbose = 1, scale = FALSE, genes_as_rows = FALSE)

m_names <- read.table("/home/tomap/Rambow_Markers_rows_Smaller.txt", sep = ' ')
x <- subset(counts, rownames(counts) %in% m_names[,1])
x <- x[order(rownames(x)),]
#x <- x[which(rowSums(x) != 0),]
#x <- x[,which(colSums(x) != 0)]

annot <- df$'1'

umap_model <- umap(t(x), random_state = 42, transform_seed = 69)

plot <- ggplot(data = data.frame(umap_model$layout), mapping = aes(x=X1,y=X2, colour = as.factor(annot))) + geom_point() + scale_colour_discrete(drop=FALSE) + labs(title = "Unsupervised Clustering with Autoencoder Preprocessing")

annot <- df$'5'
plot <- ggplot(data = data.frame(umap_model$layout), mapping = aes(x=X1,y=X2, colour = as.factor(annot))) + geom_point() + scale_colour_discrete(drop=FALSE) + labs(title = "Unsupervised Clustering with Autoencoder Preprocessing")
