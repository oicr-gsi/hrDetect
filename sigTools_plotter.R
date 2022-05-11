library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
hrd_score_file <- args[2]
sig_file <- args[3]

hrd_raw = read.table(hrd_score_file,
           sep = "\t",header = FALSE,
           stringsAsFactors = FALSE,check.names = FALSE)
hrd <- as.data.frame(t(hrd_raw$V2))
names(hrd) <- hrd_raw$V1

png(paste(sample_name,".sigtools.hrd.png",sep=""), width = 350, height = 150,type="cairo")

ggplot(hrd,
       aes(y="")) + 
  geom_vline(xintercept = 0.5,color="grey",linetype="dashed")+
  geom_point(aes(x=HRD_median), shape=4, size=6) + 
  geom_point(aes(x=HRD_point), shape=1, size=5) + 
  geom_errorbar(aes(xmin=HRD_low_quant, xmax=HRD_top_quant), width=.1) +
  theme_bw(base_size=15) + labs(y="Sample",x="HR Proficient            HR Deficient",color="VAF (<)") + 
  xlim(0,1) + guides(alpha="none")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_color_manual(values=c("#65bc45","#000000","#0099ad")) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=.5)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

sigs = read.table(sig_file,
                  sep = "\t",header = TRUE,
                  stringsAsFactors = FALSE,check.names = FALSE)
sigs$sigs <- rownames(sigs)

png(paste(sample_name,".sigtools.sigs.png",sep=""), width = 350, height = 150,type="cairo")

ggplot(sigs,aes(y="",x=sig_weight_rel_adj,fill=sigs)) + 

  geom_bar(stat="identity") + #xlim(0,1) +
  geom_text(aes(label = sigs),colour = "black",  stat="identity", position = position_stack(vjust = 0.5),vjust=-4.8) +
  labs(title="",y="") + 
  theme_bw()+ guides(fill="none") +
  xlim(0,1)+
  facet_grid(.~sig_type,scales = "free",space="free",switch = "both")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=colorRampPalette(c(rgb(101/255, 188/255, 69/255),rgb(0, 0, 0), rgb(0/255, 153/255, 173/255)), alpha = TRUE)(length(unique(sigs$sigs))))

dev.off()
  
