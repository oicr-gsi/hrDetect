library(data.table)
library(tidyverse)
library(cowplot)

setwd('/Volumes/')


test1 <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.MAFtest1.txt')
test1_dcast <- dcast(test1,formula=V1+V2+V3~V4)

#sampleSummaries <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.summaries.txt',header=F)
#sampleSummaries_dcast <- dcast(sampleSummaries,formula=V1+V2~V3)

TGL_brca2 <- fread('cgi/scratch/fbeaudry/sigTools_test/BRCA_mutations.TGL62.2.txt',header=FALSE)
PASS01_brca <- fread('cgi/scratch/fbeaudry/sigTools_test/BRCA_mutations.PASS01.txt',header=FALSE)
PASS01.2_brca <- fread('cgi/scratch/fbeaudry/sigTools_test/BRCA_mutations.PASS01.2.txt',header=FALSE)

brca <- rbind.data.frame(TGL_brca2,PASS01_brca,PASS01.2_brca) %>% select(V1,V2,V4,V5)

brca_dcast <- dcast(brca,formula=V4+V5~V1,value.var="V2")
allSamples <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.samples.txt',header=F)
brca_dcast <- left_join(allSamples,brca_dcast,by=c("V1"="V5","V2"="V4"))


test1_dcast_MAF5 <- test1_dcast %>% filter(V3 %in% c("5"))
                       
test1_dcast_MAF5$Study <-  factor(test1_dcast_MAF5$V1,  levels = c("TGL62","PASS01"))
brca_dcast$Study  <-   factor(brca_dcast$V1,  levels = c("TGL62","PASS01"))

test1_dcast_MAF5$Samples <- 
  factor(test1_dcast_MAF5$V2,  levels = test1_dcast_MAF5$V2[order(test1_dcast_MAF5$HRD_median)])
brca_dcast$Samples <- 
  factor(brca_dcast$V2, levels = test1_dcast_MAF5$V2[order(test1_dcast_MAF5$HRD_median)])

pd <- position_dodge(0.5)

HRD_plot <- 
ggplot(test1_dcast_MAF5,
       aes(x=Samples,color=Study)) + 
  geom_point(aes(y=HRD_median), shape=4, position=pd,size=3) + 
  geom_point(aes(y=HRD_point), shape=1, position=pd,size=2) + 
  
  geom_errorbar(aes(ymin=HRD_low_quant, ymax=HRD_top_quant), width=.1, position=pd) +
  theme_bw() + labs(x="Sample",y="Probability of HRD",color="Study") + 
  ylim(0,1) + guides(alpha="none")+
 # scale_shape_manual(values=c(4,1)) + 
 # facet_grid(.~as.factor(Study),scale="free",space = "free")+
   theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=c("#65bc45","#000000","#0099ad"))

BRCA_info_plot <- 
ggplot(brca_dcast ) + 
  geom_point(aes(x=Samples,y="BRCA1",color=BRCA1),size=6,shape=15) +
  geom_point(aes(x=Samples,y="BRCA2",color=BRCA2),size=6,shape=15) +
 # facet_grid(.~as.factor(Study),scale="free",space = "free")+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0)) +
  labs(x="Sample",y="",size="VAF",color="Mutation") +scale_color_discrete(na.value="white")


plot_grid(HRD_plot, BRCA_info_plot, ncol = 1, rel_heights = c(2, 1),align = 'v')




VAFvar_OCT <- fread('cgi/scratch/fbeaudry/sigTools_test/OCT_010958.VAFvar.txt')
VAFvar_PANX <- fread('cgi/scratch/fbeaudry/sigTools_test/PANX_1315.VAFvar.txt')
VAFvar <- rbind.data.frame(
  cbind.data.frame(VAFvar_OCT,"sample"="OCT_010958"),
  cbind.data.frame(VAFvar_PANX,"sample"="PANX_1315")
  
)

VAFvar_dcast <- dcast(VAFvar,formula=V3+V4+sample~V1,value.var="V2")

ggplot(VAFvar_dcast,aes(x=as.factor(V3/100),color=as.factor(V4/100))) + 
  geom_point(aes(y=HRD_median), shape=4, position=pd,size=3) + 
  geom_point(aes(y=HRD_point), shape=1, position=pd,size=2) + 

  geom_errorbar(aes(ymin=HRD_low_quant, ymax=HRD_top_quant), width=.1, position=pd) +
  theme_bw(base_size=18) + 
  facet_grid(.~as.factor(sample),scale="free",space = "free")+
  
  labs(x="INDEL VAF",y="Probability of HRD",color="SNP VAF") + 
  ylim(0,1) 

test_VAF05_summary <- test1_dcast %>% filter(V3 == 5) %>% left_join(sampleSummaries_full,by=c("V2"="V2"))
test_VAF05_summary$HRD_confidence <- test_VAF05_summary$HRD_top_quant - test_VAF05_summary$HRD_low_quant
test_VAF05_summary$BRCA <- ">1 BRCA"
test_VAF05_summary$BRCA[is.na(test_VAF05_summary$BRCA1) & is.na(test_VAF05_summary$BRCA2)] <- "No BRCA"

plot_grid(

ggplot(test_VAF05_summary,aes(x=purity,y=HRD_point,color=Study.y)) + 
  geom_point(size=2) + facet_grid(.~BRCA) +theme_bw(base_size=18) +
  labs(x="Purity",y="Probability HRD",color="Study")+
  scale_color_manual(values=c("#65bc45","#0099ad")) #+ stat_smooth(method = "lm",se=FALSE )
,
ggplot(test_VAF05_summary,aes(x=purity,y=HRD_confidence,color=Study.y)) + 
  geom_point(size=2) + facet_grid(.~BRCA) +theme_bw(base_size=18) +
  labs(x="Purity",y="HRD Confidence (size of C.I.)",color="Study")+
  scale_color_manual(values=c("#65bc45","#0099ad")) #+ stat_smooth(method = "lm",se=FALSE )
,ncol = 1, align = 'v')
purity_test <- lm(data=test_VAF05_summary,formula=purity~HRD_point*BRCA)
summary(purity_test)




ggplot(test_VAF05_summary %>% filter(V2 %in% c("OCT_010221")),
       aes(x=V2)) + 
  geom_point(aes(y=HRD_median), shape=4, position=pd,size=6) + 
  geom_point(aes(y=HRD_point), shape=1, position=pd,size=5) + 
  
  geom_errorbar(aes(ymin=HRD_low_quant, ymax=HRD_top_quant), width=.1, position=pd) +
  theme_bw(base_size=25) + labs(x="Sample",y="HR Proficient            HR Deficient",color="VAF (<)") + 
  ylim(0,1) + guides(alpha="none")+
  # scale_shape_manual(values=c(4,1)) + 
 # facet_grid(.~as.factor(Study),scale="free",space = "free")+
  theme(axis.text.x = element_text(angle = -75, vjust = 0.5, hjust=0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=c("#65bc45","#000000","#0099ad")) +
theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=.5)) 
  
