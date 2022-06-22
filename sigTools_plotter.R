library(ggplot2)
library(jsonlite)
library(cowplot)

#
args = commandArgs(trailingOnly=TRUE)
jsonIn <- args[1]

#test args
setwd('/Volumes/')
jsonInLoc  <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309.signatures.json" 

JSONin <- read_json(jsonInLoc, simplifyVector = FALSE)

sample_name <- JSONin$Sample
#####HRD-call figure####

hrd <- as.data.frame(t(c(JSONin$HRD$HRDetect$point,unlist(JSONin$HRD$HRDetect$quantiles))))
names(hrd) <- c("HRD_point","HRD_low_quant","HRD_median","HRD_top_quant")

hrd$HRD <- "HRD"

modelWeights <- as.data.frame(t(as.data.frame(JSONin$HRD$HRDetect$modelWeights)))
rownames(modelWeights)[5] <- "LOH"
modelWeights$params <- c("del.mh","SBS3","stDUPs","sDELs","LOH","SBS8")

svg(filename = paste(sample_name,".sigtools.hrd.svg",sep=""),width = 10, height = 6,type="cairo")

cowplot::plot_grid(
  
  ggplot(modelWeights,aes(x=V1,y=params)) +
    geom_vline(xintercept = 0,color="grey",linetype="dashed")+
    geom_point( shape="|", size=10) + xlim(-5,5) + theme_bw(base_size=15) +
    theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
    labs(y="") + 
    
    theme(#axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
  ,
  ggplot(hrd,
         aes(y=HRD)) + 
    geom_vline(xintercept = 0.5,color="grey",linetype="dashed")+
    geom_point(aes(x=HRD_median), shape=4, size=6) + 
    geom_point(aes(x=HRD_point), shape=1, size=5) + 
    geom_errorbar(aes(xmin=HRD_low_quant, xmax=HRD_top_quant), width=.1) +
    theme_bw(base_size=15) + 
    labs(y="",x="HR Proficient               HR Deficient") + 
    xlim(0,1) + guides(alpha="none")+
    
    scale_color_manual(values=c("#65bc45","#000000","#0099ad")) +
    theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),
  
  ncol = 1,align="tblr",rel_heights = c(2:1)
  
)

dev.off()


####signatures figure####

sigs <- as.data.frame(unlist(JSONin$SNV$classic_sigs))
names(sigs) <- "exposures"
sigs$weights <- sigs$exposures/(sum(sigs$exposures))
sigs$sig <- rownames(sigs)
sigs$SBS <- gsub(x = sigs$sig ,pattern = "Signature",replacement = "SBS.")
sigs_noZero <- sigs[sigs$exposures != 0 & sigs$sig != "unassigned",]

svg(paste(sample_name,".sigtools.sigs.svg",sep=""), width = 10, height = 6,type="cairo")

ggplot(sigs_noZero,aes(y="",x=weights,fill=SBS)) + 

  geom_bar(stat="identity") + 
  geom_text(aes(label = SBS),colour = "white",  stat="identity", position = position_stack(vjust = 0.5)) +
  labs(title="",y="") + 
  theme_bw()+ guides(fill="none") +
  xlim(0,1)+
 # facet_grid(.~sig_type,scales = "free",space="free",switch = "both")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=colorRampPalette(c(rgb(101/255, 188/255, 69/255),rgb(0, 0, 0), rgb(0/255, 153/255, 173/255)), alpha = TRUE)(length(unique(sigs_noZero$sig))))

dev.off()
  
