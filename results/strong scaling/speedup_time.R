setwd('YourPath/D-EE/results/strong scaling')
librar(ggplot2)
# strongup ratio --------------------------------------------------

tt1 <- c(1:4)
tt2 <- c(20976.7,11480.98,8501.257,7272.3)
tt3 <- c(1,tt2[1]/tt2[2],tt2[1]/tt2[3],tt2[1]/tt2[4])
tt4 <- 1:4
speedup_ratio <- data.frame(num=c(tt1, tt1), time=log2(c(tt3, c(1,2,4,8))), label=c(rep('real',4), rep('ideal',4)))

pdf('speedup ratio.pdf', height = 4, width = 4.5)
theme_set(theme_light())
ggplot(data=speedup_ratio)+
  geom_point(aes(x=num,y=time, group=label ,color = label))+
  geom_line(aes(x=num,y=time, group=label, color=label))+
  scale_x_continuous(breaks = c(1:4),labels = c("500",'1000','2000','4000'), limits = c(0.5,4.5))+
  scale_y_continuous(breaks = c(0,1,2,3), labels = c(1,2,4,8))+
  xlab("Number of Processes")+ylab('Speedup Ratio')+
  geom_label(aes(x=c(tt1,tt1),y=c(log2(tt3),log2(tt3))),label = as.character(round(c(tt3,tt3),digits = 4)))+
  theme(axis.title = element_text(size = 18, face='bold'), legend.position = c(0.9,0.15),
        legend.title = element_blank(), plot.title = element_text(hjust=0.5, face = 'bold'),
        axis.text = element_text(size = 15, face = 'bold'), 
        legend.background = element_rect(linetype = 'solid', colour = 'black', fill = 'gray95'),
        legend.text  =element_text(size = 10))+
  scale_color_manual(values = c('darkgreen','deepskyblue3'))
dev.off()


# speedup efficiency ------------------------------------------------------

speedup_effi <- data.frame(num=tt1, efficiency = tt3/c(1,2,4,8))
pdf('speedup efficiency.pdf', height = 4, width = 4.5)
theme_set(theme_light())
ggplot(data = speedup_effi, aes(x=num, y=efficiency))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = c(1:4),labels = c("500",'1000','2000','4000'), limits = c(0.5,4.5))+
  xlab("Number of Processes")+ylab('Parallel Efficiency')+
  theme(axis.title = element_text(size = 18,face='bold'),
        axis.text = element_text(size = 15, face = 'bold'))+
  geom_label(aes(x=num,y=efficiency,label = as.character(round(efficiency,digits = 4))))
dev.off()



# Time comsumption vs #Process --------------------------------------------

num_proc <- c(500,1000,2000)
N1_time <- c(5.8247e+01, 8.9404e+01,1.5533e+02)
N5_time <- c(5.9533e+02, 5.5531e+02, 5.3668e+02)
N10_time <- c(2.4125e+03, 1.6736e+03, 1.3243e+03)

ipsc_size_time <- data.frame()
N1_time <- N1_time/60
N5_time <- N5_time/60
N10_time <- N10_time/60
ipsc_size_time[1:9, 1] <- rep(c(10000,50000,100000), each =  3)
ipsc_size_time[1:9, 2] <- rep(num_proc, times = 3)
ipsc_size_time[1:9, 3] <- c(N1_time, N5_time, N10_time)
colnames(ipsc_size_time) <- c('size', 'NumProc', 'time_min')

print(ipsc_size_time)

ipsc_size_time[,1] <- rep(1:3, each=3)
ipsc_size_time$NumProc <- as.factor(ipsc_size_time$NumProc)

pdf('time_vs_size.pdf', height = 3, width = 4.3)
theme_set(theme_light())
ggplot(data = ipsc_size_time, aes(x=size, y=time_min, group=NumProc, color=NumProc))+
  geom_point()+geom_line()+
  xlab('Sample size')+ylab("Time/Min")+
  scale_color_manual(name = 'Number of\nprocesses', values = c('red3', 'deepskyblue','darkorange'))+
  scale_x_continuous(breaks = c(1:3), labels = c('10k', '50k', '100k'))+
  theme(axis.title = element_text(size = 15, face='bold'),
        axis.text = element_text(size = 13, face = 'bold'), 
        legend.background = element_rect(linetype = 'solid', colour = 'black', fill = 'gray95'),
        legend.text  =element_text(size = 10))
dev.off()







