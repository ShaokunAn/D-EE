# D-EE gene expression heatmap
all_gene <- ipsc.seurat@assays$RNA@data@Dimnames[[1]]
gene_name <- c('Sox2', 'Sox4', 'Nanog')
for(i in 1:length(gene_name)){
  gene_id <- which(all_gene==gene_name[i])
  plot_data <- data.frame(dee1 = ipsc.dee.2d[,1], dee2 = ipsc.dee.2d[,2],
                          gene = ipsc.seurat@assays$RNA@data[gene_id,])
  plot_data2 <- plot_data[order(plot_data$gene),]
  pdf(paste('ipsc_dee_', gene_name[i], '.pdf', sep = ''), height = 6, width = 7)
  theme_set(theme_bw())
  ggplot(data = plot_data2, aes(x=dee1, y=dee2)) + geom_point(aes(color = gene), size = 0.3)+
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 15), 
          legend.title = element_text(size = 15), legend.text = element_text(size = 13))+
    labs(x = 'DEE1', y='DEE2')+
    scale_color_gradient2(low = 'grey',mid = 'red', high = 'black', midpoint = mean(range(plot_data2$gene)), name = gene_name[i])
  dev.off()
}

# D-TSEE gene expression heatmap
for(i in 1:length(gene_name)){
  gene_id <- which(all_gene==gene_name[i])
  plot_data <- data.frame(dtsee1 = ipsc.dtsee.2d[,1], dtsee2 = ipsc.dtsee.2d[,2],
                          gene = ipsc.seurat@assays$RNA@data[gene_id,])
  plot_data2 <- plot_data[order(plot_data$gene),]
  pdf(paste('ipsc_dtsee_', gene_name[i], '.pdf', sep = ''), height = 6, width = 7)
  theme_set(theme_bw())
  ggplot(data = plot_data2, aes(x=dtsee1, y=dtsee2)) + geom_point(aes(color = gene), size = 0.3)+
    theme(axis.title = element_text(size = 20, colour = 'black'), axis.text = element_text(size = 15), 
          legend.title = element_text(size = 15), legend.text = element_text(size = 13))+
    labs(x = 'DTSEE1', y='DTSEE2')+
    scale_color_gradient2(low = 'grey',mid = 'red', high = 'black', midpoint = mean(range(plot_data2$gene)), name = gene_name[i])
  dev.off()
}
