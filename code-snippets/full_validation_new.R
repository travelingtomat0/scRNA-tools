# This program conducts the cell assignments for each subdirectory

# To keep the command running, according to
# https://medium.com/@arnab.k/how-to-keep-processes-running-after-ending-ssh-session-c836010b26a3
# use ssh, then screen, then can disconnect with ctrl-A & ctrl-D and later can resume session with screen -r

# TODO: Impute by -1, if marker gene was not counted in analysis?
# TODO: MAHEFOG (5)

library(hdf5r)
library(rhdf5)
library(ggplot2)
library(monocle)
library(garnett)
library(org.Hs.eg.db)

# Build some basic functions used for validation.

# something wrong in this function with length of tab...

meanCell <- function(tab) {
 temp <- integer(nrow(tab))
 for(i in 1:ncol(tab)) {
  temp <- temp + tab[,i]
 }
 return(round((1/ncol(tab)) * temp))
}

# cts is for instance counts[,invasive], i.e. type-spec. count data.
# gives avg. expression 
intraclass_sim <- function(cts, avg_cell) {
	res <- 0
	for(i in 1:ncol(cts)) {
		res = res + cor(avg_cell, cts[,i], method = "pearson")
	}
	return((1/ncol(cts))*res)
}

# Build some base global variables
# The assumption is that this file is run from the "~/infercnv_annotations" folder
# Which contains all relevant directories and files in the resp. sub-directories.
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
m_names <- read.table("/home/tomap/Rambow_Markers_rows_Smaller.txt", sep = ' ')

avg_invasive_set <- data.frame(integer(nrow(m_names)))
rownames(avg_invasive_set) <- m_names[,1]
avg_pigmented_set <- data.frame(integer(nrow(m_names)))
rownames(avg_pigmented_set) <- m_names[,1]
avg_ncsc_set <- data.frame(integer(nrow(m_names)))
rownames(avg_ncsc_set) <- m_names[,1]
avg_smc_set <- data.frame(integer(nrow(m_names)))
rownames(avg_smc_set) <- m_names[,1]

information_vec <- c()
information_names <- c()
total_inv = 0
total_pig = 0
total = 0

# Iteration over all tumors in the directory
for(i in 1:length(subdirs)) {
# Garbage collection after last iteration
gc()

# Check, if sample was actually analysed.
if(!file.exists(paste(subdirs[i],"/warning.txt",sep=""))) {

print(paste("Begin work in directory", subdirs[i], sep = " "))

exc <- try(dir.create(paste(subdirs[i],"/validation", sep = "")))

# read in counts data and previous annotations.
infile <- H5Fopen(paste(subdirs[i],"/count_data.h5", sep = ""))
counts <- infile$raw_counts
rownames(counts) <- make.unique(infile$gene_attrs$gene_names)
colnames(counts) <- infile$cell_attrs[[1]]


print("Finished reading in files")

x <- subset(counts, rownames(counts) %in% m_names[,1])
x <- x[order(rownames(x)),]

ass <- data.frame(c())

FLAG1 <- FALSE
if(file.exists(paste(subdirs[i],"/results/assignment3_cell_types.txt", sep = ""))) {
	ass <- try(read.table(paste(subdirs[i],"/results/assignment3_cell_types.txt", sep = ""), sep= " "))
	if(inherits(ass, "try-error")) {
		print(paste("Directory",subdirs[i],"doesn't contain assignments", sep = " "))
	} else {
		FLAG1 <- TRUE
	}
}

# TEST if assignments exist. If not, then dont get past the flag.
if(FLAG1 && dir.exists(paste(subdirs[i],"/results", sep = ""))) {
names <- data.frame(rownames(ass))

invasive <- names[ass$Cell.Type == "Invasive",]
ncsc <- names[ass$Cell.Type == "NCSC",]
pigmented <- names[ass$Cell.Type == "Pigmented",]
smc <- names[ass$Cell.Type == "SMC",]
# unassigned <- try(x[,ass$Cell.Type == "Unassigned"])

counts <- counts[,order(colnames(counts))]
c <- new_cell_data_set(counts)
#c <- estimateSizeFactors(c)
marker_check <- check_markers(c, "~/new_infercnv_annotations/garnett_markup.txt", db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
pdf(paste(subdirs[i],"/results","/marker_file.pdf", sep = ""))
plot_markers(marker_check)
dev.off()
write.table(marker_check, file = paste(subdirs[i],"/results","/marker_information.txt", sep = ""), sep = ',')

# This saves the marker check in the results section.
# To save a file that checks marker gene expression for the specific clusters, use the commented code below.
# check_markers(new_cell_data_set(data.matrix(counts[,invasive])), "~/new_infercnv_annotations/garnett_markup.txt", db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
# TODO: Can save a table containing ambiguity, marker score etc. Its all in the marker_check object.
# I.e. use code: write.table(data = marker_check, file = paste(subdirs[i],"/results","/marker_information.txt", sep = ""), sep = ',')

# Write information-vector for Agnieszka:
print("Write information for Agnieszka...")
information_vec <- c(information_vec, length(invasive)/(length(invasive) + length(pigmented)), length(invasive)/nrow(ass), length(pigmented)/nrow(ass))
#information_vetcor <- c(information_vector, nrow(invasive)/nrow(annotations))
information_names <- c(information_names, paste(substring(subdirs[i],3,9), "Invasive_Pigmented", sep = "_"), paste(substring(subdirs[i],3,9), "Invasive_All", sep = "_"), paste(substring(subdirs[i],3,9), "Pigmented_All", sep = "_"))
total_inv <- total_inv + length(invasive)
total_pig <- total_pig + length(pigmented)
total <- total + length(invasive) + length(pigmented) + length(ncsc) + length(smc)

###########################################################################################
# Actual Validation (intra-sample)

pig_empty <- FALSE
inv_empty <- FALSE
ncsc_empty <- FALSE
smc_empty <- FALSE

print("Calculate avg cells...")
if(length(invasive) != 0) {
	if(length(invasive) == 1) {
		avg_invasive_cell <- counts[,invasive]
		inv_empty <- TRUE
	} else
		avg_invasive_cell <- meanCell(counts[, invasive])
		pdf(paste(subdirs[i],"/results","/marker_file_invasive.pdf", sep = ""))
		c <- new_cell_data_set(data.matrix(counts[,invasive]))
		marker_check <- check_markers(c, "~/new_infercnv_annotations/garnett_markup.txt", db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
		plot_markers(marker_check)
		dev.off()
		write.table(marker_check, file = paste(subdirs[i],"/results","/marker_invasive_information.txt", sep = ""), sep = ',')
	}
} else {
	avg_invasive_cell <- data.frame(c(-2*integer(nrow(counts))))
	rownames(avg_invasive_cell) <- rownames(counts)
	inv_empty <- TRUE
}
if(length(ncsc) != 0) {
	if(length(ncsc) == 1) {
		avg_ncsc_cell <- counts[,ncsc]
		ncsc_empty <- TRUE
	} else {
		avg_ncsc_cell <- meanCell(counts[,ncsc])
		pdf(paste(subdirs[i],"/results","/marker_file_ncsc.pdf", sep = ""))
		c <- new_cell_data_set(data.matrix(counts[,ncsc]))
		marker_check <- check_markers(c, "~/new_infercnv_annotations/garnett_markup.txt", db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
		plot_markers(marker_check)
		dev.off()
		write.table(marker_check, file = paste(subdirs[i],"/results","/marker_ncsc_information.txt", sep = ""), sep = ',')
	}
} else {
	avg_ncsc_cell <- data.frame(c(-2*integer(nrow(counts))))
	rownames(avg_ncsc_cell) <- rownames(counts)
	ncsc_empty <- TRUE
}
if(length(pigmented) != 0) {
	avg_pigmented_cell <- meanCell(counts[,pigmented])
	pdf(paste(subdirs[i],"/results","/marker_file_pigmented.pdf", sep = ""))
	c <- new_cell_data_set(data.matrix(counts[,pigmented]))
	marker_check <- check_markers(c, "~/new_infercnv_annotations/garnett_markup.txt", db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
	plot_markers(marker_check)
	dev.off()
	write.table(marker_check, file = paste(subdirs[i],"/results","/marker_pigmented_information.txt", sep = ""), sep = ',')

} else {
	avg_pigmented_cell <- data.frame(c(-2*integer(nrow(counts))))
	rownames(avg_pigmented_cell) <- rownames(counts)
	pig_empty <- TRUE
}
if(length(smc) != 0) {
	if(length(smc) == 1) {
		avg_smc_cell <- counts[,smc]
		smc_empty <- TRUE
	} else {
		avg_smc_cell <- meanCell(counts[,smc])
		pdf(paste(subdirs[i],"/results","/marker_file_smc.pdf", sep = ""))
		c <- new_cell_data_set(data.matrix(counts[,smc]))
		marker_check <- check_markers(c, "~/new_infercnv_annotations/garnett_markup.txt", db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")
		plot_markers(marker_check)
		dev.off()
		write.table(marker_check, file = paste(subdirs[i],"/results","/marker_smc_information.txt", sep = ""), sep = ',')
	}
} else {
	avg_smc_cell <- data.frame(c(-2*integer(nrow(counts))))
	rownames(avg_smc_cell) <- rownames(counts)
	smc_empty <- TRUE
}
# avg_unassigned_cell <- meanCell(unassigned)

print("Finished claculating avg cells.")

df <- data.frame(c(data.frame(avg_invasive_cell),data.frame(avg_pigmented_cell),data.frame(avg_ncsc_cell), data.frame(avg_smc_cell)))
rownames(df) <- rownames(counts)
colnames(df) <- c("Avg_Invasive", "Avg_Pigmented", "Avg_NCSC", "Avg_SMC")
df <- df[order(rownames(df)),]

print(subdirs)

write.table(df, file = paste(subdirs[i],"/validation/avg_cell_fullgene.txt", sep = ""), sep = ",")

mf <- data.frame(df[rownames(df) %in% m_names[,1],])
mf <- mf[order(rownames(mf)),]
write.table(mf, file = paste(subdirs[i],"/validation/avg_cell_markergene.txt", sep = ""), sep = ",")

# TODO: Impute by -1, if marker gene was not counted in analysis?
# Could use function impute(avg_invasive_cell, m_names[,1])
# TODO: Important question.. Use normalized values? Then we would need to normalize the whole counts
# matrix before computing and comparing average cells.
# TODO: Also: normalization?
#print("Calculate avg.-sets...")

#avg_invasive_set <- data.frame(c(avg_invasive_set, avg_invasive_cell))
#avg_ncsc_set <- c(avg_ncsc_set, avg_ncsc_cell)
#avg_pigmented_set <- c(avg_pigmented_set, avg_pigmented_cell)
#avg_smc_set <- c(avg_smc_set, avg_smc_cell)

#TODO: Top is for full set of genes. Bottom is for reduced to marker genes.
#avg_invasive_set <- c(avg_invasive_set, mf$Avg_Invasive)
#avg_ncsc_set <- c(avg_ncsc_set, mf$Avg_NCSC)
#avg_pigmented_set <- c(avg_pigmented_set, mf$Avg_Pigmented)
#avg_smc_set <- c(avg_smc_set, mf$Avg_SMC)

# Claculate average intraclass correlation (only withmarker genes for now)
# ONLY one expressed cell counts as not expressed!
inv_corr <- ifelse(inv_empty,NA, intraclass_sim(x[,invasive],mf$Avg_Invasive))
pig_corr <- ifelse(pig_empty, NA,intraclass_sim(x[,pigmented],mf$Avg_Pigmented))
ncsc_corr <- ifelse(ncsc_empty, NA,intraclass_sim(x[,ncsc],mf$Avg_NCSC))
smc_corr <- ifelse(smc_empty, NA,intraclass_sim(x[,smc],mf$Avg_SMC))

write.table(c(paste("Invasive correlation", inv_corr, sep = ": "), paste("Pigmented correlation", pig_corr, sep = ": "), paste("NCSC correlation", ncsc_corr, sep = ": "), paste("SMC Correlation", smc_corr, sep = ": "), pig_empty, inv_empty, ncsc_empty, smc_empty), file = paste(subdirs[i],"/validation/",substring(subdirs[i],3,10),"intraclass_correlations.txt", sep = ""), sep = ",")

inv_piga <- ifelse(inv_empty,NA, intraclass_sim(x[,invasive],mf$Avg_Pigmented))
inv_smca <- ifelse(inv_empty,NA, intraclass_sim(x[,invasive],mf$Avg_SMC))
pig_inva <- ifelse(pig_empty, NA,intraclass_sim(x[,pigmented],mf$Avg_Invasive))
pig_smca <- ifelse(pig_empty, NA,intraclass_sim(x[,pigmented],mf$Avg_SMC))
smc_inva <- ifelse(smc_empty, NA,intraclass_sim(x[,smc],mf$Avg_Invasive))
smc_piga <- ifelse(smc_empty, NA,intraclass_sim(x[,smc],mf$Avg_Pigmented))

write.table(c(c("Inv-PigAvg", "Inv-SmcAvg", "Pig-InvAvg", "Pig-SmcAvg", "Smc-InvAvg", "Smc-PigAvg"), c(inv_piga, inv_smca, pig_inva, pig_smca, smc_inva, smc_piga)), file = paste(subdirs[i],"/validation/",substring(subdirs[i],3,10),"interclass_correlations.txt", sep = ""), sep = ",")
}
}

information_vec <- c(information_vec, total_inv/total_pig, total_pig/total, total_inv/total)
information_names <- c(information_names, "Invasive-Pigmented (Total)", "Pigmented-Total", "Invasive-Total")
information_vector <- data.frame(information_vec)
rownames(information_vector) <- information_names
colnames(information_vector) <- c("Fraction")
write.table(information_vector, "information_for_agnieszka.txt", sep = '\t')

warnings()

