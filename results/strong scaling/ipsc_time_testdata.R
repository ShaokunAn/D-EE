# generate test data to obtain the time consumption of various numbers of data  
setwd('YourPath/D-EE/results/strong scaling')
ipsc_afterpca_all <- read.csv('ipsc_afterpca.csv', header = FALSE)
ipsc_afterpca_1w <- ipsc_afterpca_all[1:10000, ]
ipsc_afterpca_5w <- ipsc_afterpca_all[1:50000,]
ipsc_afterpca_10w <- ipsc_afterpca_all[1:100000,]
write.table(ipsc_afterpca_1w, file = 'ipsc_afterpca_1w.csv', 
            col.names = FALSE, row.names = FALSE, sep = ',')
write.table(ipsc_afterpca_5w, file = 'ipsc_afterpca_5w.csv', 
            col.names = FALSE, row.names = FALSE, sep = ',')
write.table(ipsc_afterpca_10w, file = 'ipsc_afterpca_10w.csv', 
            col.names = FALSE, row.names = FALSE, sep = ',')
