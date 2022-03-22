library(data.table)
library(tidyverse)

setwd('/Volumes/')
pd <- position_dodge(0.25)


test1 <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.MAFtest1.txt')
test1_dcast <- dcast(test1,formula=V1+V2+V3~V4)

sampleSummaries <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.summaries.txt',header=F)
sampleSummaries_dcast <- dcast(sampleSummaries,formula=V1+V2~V3)

test1_dcast$Samples <- 
  factor(test1_dcast$V2, levels = sampleSummaries_dcast$V2[order(sampleSummaries_dcast$purity)])

HRD_plot <- 
  
ggplot(test1_dcast,aes(x=Samples,color=as.factor(V3/100))) + 
  geom_point(aes(y=HRD_median), shape=4, position=pd) + 
  geom_point(aes(y=HRD_point), shape=1, position=pd) + 
  
  geom_errorbar(aes(ymin=HRD_low_quant, ymax=HRD_top_quant), width=.1, position=pd) +
  theme_bw() + labs(x="Sample",y="Probability of HRD",color="VAF (<)") + 
  ylim(0,1) + 
 # scale_shape_manual(values=c(4,1)) + 
  facet_grid(.~as.factor(V1),scale="free",space = "free")+
   theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

sampleSummaries_dcast$Samples <- 
  factor(sampleSummaries_dcast$V2, levels = sampleSummaries_dcast$V2[order(sampleSummaries_dcast$purity)])

BRCA_info_plot <- 
ggplot(sampleSummaries_dcast %>% filter(V2 != "PANX_1312")) + 
  geom_point(aes(x=Samples,y="BRCA1",shape=BRCA1_mutation),size=5) +
geom_point(aes(x=Samples,y="BRCA2",shape=BRCA2_mutation),size=5) +
  facet_grid(.~as.factor(V1),scale="free",space = "free")+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0)) +
  labs(x="Sample",y="",size="VAF",shape="Mutation") 

library(cowplot)

plot_grid(HRD_plot, BRCA_info_plot, ncol = 1, rel_heights = c(2, 1),align = 'v')
