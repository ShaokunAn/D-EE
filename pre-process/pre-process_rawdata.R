# read raw datasets and pre-process them
# notice the directory paths of the raw datasets need to be set correctly and manually

install.packages('R.matlab')
install.packages('ggplot2')
install.packages('Seurat')
install.packages('shape')
library(R.matlab)
library(ggplot2)
library(Seurat)
library(shape)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('rhdf5')
library(rhdf5)

DIR <- 'YourPath/D-EE/'

# pre-process of PHATE dataset --------------------------------------------

phate_raw <- readMat(paste(DIR, 'pre-process/PHATE_Data.mat', sep=''))
phate_raw <- phate_raw$M
after_pca <- function(data){
  data_guiyi <- (data-min(data))/(max(data)-min(data))
  data_cov <- cov(data_guiyi)
  data_eigen <- eigen(data_cov)
  lambdanum <- min(which(abs(diff(data_eigen$values,2))<1e-3))
  data_pca <- data_guiyi%*%data_eigen$vectors[,1:lambdanum]
  return(data_pca)
}
phate_pca <- after_pca(phate_raw)
write.table(phate_pca, file = 'phate_pca.csv', row.names=FALSE, col.names = FALSE, sep = ',')


# Pre-process HSPCs dataset -----------------------------------------------

# first download the dataset from https://github.com/dpeerlab/wishbone 
# and select data in ./wishbone/data/sample_scseq_data.csv 
hspcs <- read.csv('YourPath/wishbone/data/sample_scseq_data.csv')
cell_name <- hspcs[,1]
cell_name <- as.character(cell_name)
hspcs <- hspcs[,2:2313]
hspcs <- t(hspcs)
colnames(hspcs) <- cell_name
hspcs <- CreateSeuratObject(hspcs)
hspcs <- NormalizeData(hspcs, normalization.method = 'LogNormalize')
hspcs <- FindVariableFeatures(hspcs, selection.method = 'vst', nfeatures = 2000)
hspcs <- ScaleData(hspcs)
hspcs <- RunPCA(hspcs, features = VariableFeatures(object = hspcs))
hspcs_pca <- hspcs@reductions$pca@cell.embeddings
write.table(hspcs_pca, file = 'hspcs_pca.csv', row.names=FALSE, col.names = FALSE, sep = ',')


# pre-process iPSCs data --------------------------------------------------

# first download GSE115943_RAW data and unzip it
# in command line, conduct "ls > name.txt"  then perform following codes in R
filename <- read.table("name.txt", stringsAsFactors = FALSE)

ipsc.data <- list()

for(i in 1:nrow(filename)){
  ipsc.data[[i]] <- h5read(file = filename[i,1], name = "mm10")  
}

ipsc.all.dgtmat <- list()
for(i in 1:length(ipsc.data)){
  if(length(ipsc.data[[i]]$indptr)!=length(unique(ipsc.data[[i]]$indptr))){print(i)}
  tt1 <- c()
  for(j in 1:(length(ipsc.data[[i]]$indptr)-1)){
    tt <- rep(j, each = ipsc.data[[i]]$indptr[j+1]-ipsc.data[[i]]$indptr[j])
    tt1 <- c(tt1, tt)
  }
  tt2 <- sparseMatrix(i = ipsc.data[[i]]$indices+1, j=tt1,
                      x = ipsc.data[[i]]$data, dims = c(ipsc.data[[i]]$shape[1],ipsc.data[[i]]$shape[2]))
  # sparseMatrix's index is 1-based index, not 0-based index
  tt2 <- as(tt2, 'dgTMatrix')
  tt2@Dimnames[[1]] <- ipsc.data[[i]]$gene_names
  tt2@Dimnames[[2]] <- ipsc.data[[i]]$barcodes
  # dgTMatrix's index is 0-based
  if(i==1){
    ipsc.all.dgtmat <- tt2
  }else{
    ipsc.all.dgtmat <- cbind(ipsc.all.dgtmat, tt2)
  }
  
}
rm(ipsc.data)
ipsc.all.dgtmat@Dimnames[[2]] <- paste(1:ipsc.all.dgtmat@Dim[2], ipsc.all.dgtmat@Dimnames[[2]], sep = '_')
ipsc.seurat <- CreateSeuratObject(counts=ipsc.all.dgtmat, min.cells = 50, min.features = 200)
ipsc.seurat <- NormalizeData(ipsc.seurat, normalization.method = 'LogNormalize')
ipsc.seurat <- FindVariableFeatures(ipsc.seurat, selection.method = 'vst', nfeatures = 2000)
ipsc.seurat <- ScaleData(ipsc.seurat)
ipsc.seurat <- RunPCA(ipsc.seurat, features = VariableFeatures(object = ipsc.seurat))
ipsc.pca <- ipsc.seurat@reductions$pca@cell.embeddings
colnames(ipsc.pca) <- NULL
rownames(ipsc.pca) <- NULL
write.matrix(ipsc.pca, file = 'ipsc_afterpca.csv', sep = ',')
# you can save the preprocessed data in case
# base::save.image(file = 'D-EE_data.RData')

# settle the time point of each sample in iPSCs for D-TSEE software
filtered_cellid <- ipsc.seurat@assays$RNA@counts@Dimnames[[2]]
origin_cellid <- ipsc.all.dgtmat@Dimnames[[2]]
missingid <- c()
j <- 1
for(i in 1:length(origin_cellid)){
  if(origin_cellid[i]==filtered_cellid[j]){
    j <- j+1
  }else{
    missingid <- c(missingid, i)
  }
}

# settle the time point of each sample in iPSCs for D-TSEE software
tt1 <- rep(seq(from = 0, to = 8, by = 0.5), each=2)
tt1 <- c(tt1, rep(seq(from=8.25, to=9,by=0.25), each=4))
tt1 <- c(tt1, rep(seq(from=9.5, to=18, by=0.5), each=4))
tt1 <- c(tt1, rep(100, each=4)) # two datasets have no time labels, and we use 100 to label them virtually
names(tt1) <- NULL # tt1 is the 39 various time points 


tt2 <- c() # tt2 is the number of samples in each time point
for(i in 1:length(ipsc.data)){
  tt2 <- c(tt2, ipsc.data[[i]]$shape[2])
}
cell_time <- c()
for(i in 1:length(tt2)){
  cell_time <- c(cell_time, rep(tt1[i], each=tt2[i]))
}
cell_time <- cell_time[-missingid]

cell_time_tsee <- cell_time
cell_time_tsee[which(cell_time_tsee==100)] <- 20
write.table(cell_time_tsee, file = 'ipsc_time_tsee.csv', row.names = FALSE, col.names = FALSE, sep = ',')

