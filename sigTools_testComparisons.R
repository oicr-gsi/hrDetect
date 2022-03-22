library(data.table)
library(tidyverse)

setwd('/Volumes/')

test1 <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.test1.text')

test1_dcast <- dcast(test1,formula=V1+V2+V3~V4)

test1_HRDscore <- test1 %>% filter(V4 %in% c("HRD_point","HRD_median"))
test1_HRDscore$V5 <- as.numeric(test1_HRDscore$V5)

#test1_HRDscore$Sample <- factor(test1_HRDscore$V2[order(test1_HRDscore$V5)], levels = test1_HRDscore$V2[order(test1_HRDscore$V5)])

ggplot(test1_HRDscore,aes(x=V2,y=V5,color=as.factor(V3),shape=V4)) + 
  geom_point() + theme_bw() + labs(x="Sample",y="Probability of HRD",shape="HRD") + 
  ylim(0,1) + scale_shape_manual(values=c(4,1)) + facet_grid(as.factor(V3)~.)+
   theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0))


ggplot(test1_dcast,aes(x=as.numeric(HRD_point),y=as.numeric(purity))) + geom_point()
