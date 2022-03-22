library(data.table)
library(tidyverse)

setwd('/Volumes/')

test1 <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.MAFtest1.txt')

test1_dcast <- dcast(test1,formula=V1+V2+V3~V4)

pd <- position_dodge(0.2)
ggplot(test1_dcast,aes(x=V2,color=as.factor(V3/100))) + 
  geom_point(aes(y=HRD_median), shape=4, position=pd) + 
  geom_point(aes(y=HRD_point), shape=1, position=pd) + 
  
  geom_errorbar(aes(ymin=HRD_low_quant, ymax=HRD_top_quant), width=.1, position=pd) +
  theme_bw() + labs(x="Sample",y="Probability of HRD",color="VAF") + 
  ylim(0,1) + 
 # scale_shape_manual(values=c(4,1)) + 
  facet_grid(.~as.factor(V1),scale="free",space = "free")+
   theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0))


sampleSummaries <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.summaries.txt',header=F)
sampleSummaries <- dcast(sampleSummaries,formula=V1+V2~V3)
