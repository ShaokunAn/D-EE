
# plot PHATE 2d embeddings obtained by D-EE and EE, respectively ----------
# perform D-EE in C and read it result
library(ggplot2)
setwd('YourPath/D-EE/results/consistency validation')

phate_pca_dee <- read.csv('phate_matlabX0_dee.csv', header = FALSE)
phate_pca_dee <- as.data.frame(phate_pca_dee)
phate_pca_ee <- read.csv('phate_matlabX0_ee.csv', header = FALSE)
phate_pca_ee <- as.data.frame(phate_pca_ee)
colnames(phate_pca_dee) <- c('DEE1', 'DEE2')
colnames(phate_pca_ee) <- c('EE1', 'EE2')
pdf('PHATE_DEE.pdf', height = 4, width = 5)
theme_set(theme_light())
ggplot(data = phate_pca_dee)+geom_point(aes(x=DEE1, y=DEE2), size = 0.4)+
  ggtitle("D-EE")+
  theme(axis.title = element_text(size = 18), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 13))
dev.off()

pdf('PHATE_EE.pdf', height = 4, width = 5)
theme_set(theme_light())
ggplot(data = phate_pca_ee)+geom_point(aes(x=EE1, y=EE2), size = 0.4)+
  ggtitle("EE")+
  theme(axis.title = element_text(size = 18), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 13))
dev.off()


# plot HSPCs 2d embeddings obtained by D-EE and EE, respectively ----------

hspcs_pca_dee <- read.csv('hspcs_matlabX0_dee.csv', header = FALSE)
hspcs_pca_ee <- as.data.frame(hspcs_pca_dee)
colnames(hspcs_pca_dee) <- c('DEE1','DEE2')
hspcs_pca_ee <- read.csv('hspcs_matlabX0_ee.csv', header = FALSE)
hspcs_pca_ee <- as.data.frame(hspcs_pca_ee)
colnames(hspcs_pca_ee) <- c('EE1','EE2')
pdf('HSPCs_DEE.pdf', height = 4, width = 5)
theme_set(theme_light())
ggplot(data = hspcs_pca_dee)+geom_point(aes(x=DEE1, y=DEE2), size = 0.4)+
  ggtitle("D-EE")+
  theme(axis.title = element_text(size = 18), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 13))
dev.off()

pdf('HSPCs_EE.pdf', height = 4, width = 5)
theme_set(theme_light())
ggplot(data = hspcs_pca_ee)+geom_point(aes(x=EE1, y=EE2), size = 0.4)+
  ggtitle("EE")+
  theme(axis.title = element_text(size = 18), 
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text = element_text(size = 13))
dev.off()


