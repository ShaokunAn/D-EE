library(ggplot2)
library(shape)
setwd('YourPath/D-EE/results/t-SNE_UMAP_DEE_DTSEE comparison')

# D-EE result

ipsc.dee.2d <- read.csv('ipsc_pca_dee.csv', header = FALSE)
cell_time <- read.csv('ipsc_time_tsee.csv', header = FALSE)
cell_time <- cell_time$V1

col <- intpalette(c('black','blue','cyan','green','yellow','orange','red'), 101)
mycol <- round(cell_time/unique(cell_time)[length(unique(cell_time))-1]*100+1)
mycol <- mycol[which(mycol<102)]
mycol <- col[mycol]


pdf('ipsc_dee_time.pdf', height = 6, width = 7)
# png('ipsc_dee_time.png')
par(mgp = c(3,1,0), mar = c(5,5,4,5.5)+0.1)
plot(x=ipsc.dee.2d[which(cell_time<20),1], y = ipsc.dee.2d[which(cell_time<20), 2], col = mycol,
     pch = 19, cex = 0.1, xlab = 'DEE1', ylab = 'DEE2', cex.lab = 1.3)
points(x = ipsc.dee.2d[which(cell_time==20), 1], y = ipsc.dee.2d[which(cell_time==20), 2], 
       col = 'red4', pch = 19, cex = 0.1)
colorlegend(col = col, zlim = c(0, 18), dz =3, digit = 0, posx = c(0.88, 0.91), posy = c(0.25, 0.75))
legend(x = 1.1, y = 0.85, legend = 'iPSCs', col = 'red4', bty = 'n', pch = 19)
dev.off()


# t-SNE result
# If the Seurat object of iPSCs has been removed, 
# you need to load them which has been saved in pre-process step
# or repeat the pre-process step.

# when using the FIt-SNE method, you need to compile it in advance
# for details, please refer to the instructions of FIt-SNE in https://github.com/KlugerLab/FIt-SNE
# Besides, we also provide our FIt-SNE output of iPSCs data, you can use it
# directly by the following command:
# ipsc.tsne.2d <- read.csv('ipsc_pca_fitsne.csv', header = FALSE)

ipsc.seurat <- RunTSNE(ipsc.seurat, dims = 1:50, tsne.method = 'FIt-SNE', fast_tsne_path = 'YourDir/fast_tsne', learning_rate = dim(ipsc.seurat)[2]/12, initializatinon='pca')
ipsc.tsne.2d <- ipsc.seurat@reductions$tsne@cell.embeddings
pdf('ipsc_tsne_time.pdf', height = 6, width = 7)
# png('ipsc_dee_time.png')
par(mgp = c(3,1,0), mar = c(5,5,4,5.5)+0.1)
plot(x=ipsc.tsne.2d[which(cell_time<20),1], y = ipsc.tsne.2d[which(cell_time<20), 2], col = mycol,
     pch = 19, cex = 0.1, xlab = 'tSNE1', ylab = 'tSNE2', cex.lab = 1.3)
points(x = ipsc.tsne.2d[which(cell_time==20), 1], y = ipsc.tsne.2d[which(cell_time==20), 2], 
       col = 'red4', pch = 19, cex = 0.1)
colorlegend(col = col, zlim = c(0, 18), dz =3, digit = 0, posx = c(0.88, 0.91), posy = c(0.25, 0.75))
legend(x = 1.1, y = 0.85, legend = 'iPSCs', col = 'red4', bty = 'n', pch = 19)
dev.off()


# UMAP result
# the ipsc.seurat object is still needed for running UMAP
# and we also provide our UMAP output of iPSCs data which can be used
# directly with the following command:
# ipsc.umap.2d <- read.csv('ipsc_pca_umapnb100.csv', header = FALSE)
ipsc.seurat <- RunUMAP(ipsc.seurat, dims = 1:50, n.neighbors = 100L)
ipsc.umap.2d <- ipsc.seurat@reductions$umap@cell.embeddings
pdf('ipsc_umap_time.pdf', height = 6, width = 7)
# png('ipsc_dee_time.png')
par(mgp = c(3,1,0), mar = c(5,5,4,5.5)+0.1)
plot(x=ipsc.umap.2d[which(cell_time<20),1], y = ipsc.umap.2d[which(cell_time<20), 2], col = mycol,
     pch = 19, cex = 0.1, xlab = 'UMAP1', ylab = 'UMAP2', cex.lab = 1.3)
points(x = ipsc.umap.2d[which(cell_time==20), 1], y = ipsc.umap.2d[which(cell_time==20), 2], 
       col = 'red4', pch = 19, cex = 0.1)
colorlegend(col = col, zlim = c(0, 18), dz =3, digit = 0, posx = c(0.88, 0.91), posy = c(0.25, 0.75))
legend(x = 1.1, y = 0.85, legend = 'iPSCs', col = 'red4', bty = 'n', pch = 19)
dev.off()


#D-TSEE result
ipsc.dtsee.2d <- read.csv('ipsc_afterpca_dtsee.csv', header = FALSE)
pdf('ipsc_dtsee_time.pdf', height = 6, width = 7)
# png('ipsc_dee_time.png')
par(mgp = c(3,1,0), mar = c(5,5,4,5.5)+0.1)
plot(x=ipsc.dtsee.2d[which(cell_time<20),1], y = ipsc.dtsee.2d[which(cell_time<20), 2], col = mycol,
     pch = 19, cex = 0.1, xlab = 'DTSEE1', ylab = 'DTSEE2', cex.lab = 1.3)
points(x = ipsc.dtsee.2d[which(cell_time==20), 1], y = ipsc.dtsee.2d[which(cell_time==20), 2], 
       col = 'red4', pch = 19, cex = 0.1)
colorlegend(col = col, zlim = c(0, 18), dz =3, digit = 0, posx = c(0.88, 0.91), posy = c(0.25, 0.75))
legend(x = 1.1, y = 0.85, legend = 'iPSCs', col = 'red4', bty = 'n', pch = 19)
dev.off()
