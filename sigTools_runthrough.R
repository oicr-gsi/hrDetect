##version 1.3

####packages####
#install_github('Nik-Zainal-Group/signature.tools.lib',ref='v2.1.2')
library(signature.tools.lib)
library(optparse)

library(CHORD)
library(BSgenome.Hsapiens.UCSC.hg38)
library(jsonlite)

####functions####
'%ni%' <- function(x,y)!('%in%'(x,y))

sigTools_formatter <- function(input,sampleName){
  cat_list <- list()
  names(input)[1] <- "catalogue"
  colnames(input$catalogue) <- sampleName
  cat_list[[1]] <- input$catalogue
  catalogues <- do.call(cbind,cat_list)
  return(catalogues)
}

reformat_toJSON <- function(inCol){
  catalog <-   as.data.frame(t(inCol[,1])) 
  names(catalog) <- rownames(inCol) 
  return(catalog)
}

setColNames <- function (object = nm, nm) {
  colnames(object) <- nm
  object
}

#extend_segments is not currently in use, empirically found it did not help with HRD
extend_segments <-  function(ascat.data){
  #extend disrupted segments by assuming homozygosity within missing chunk 
  #is same on both sides of missing chunk; increases # of LOH segs
  
  df_total = data.frame()
  
  for(chrom in unique(ascat.data$seg_no)){
    
    first.seg <- min(unique(ascat.data$seg_no[ascat.data$Chromosome == chrom]))
    
    start <- ascat.data$chromStart[ascat.data$seg_no == first.seg]
    end   <- ascat.data$chromEnd[ascat.data$seg_no == first.seg]
    nT    <- ascat.data$total.copy.number.inTumour[ascat.data$seg_no == first.seg]
    nB    <- ascat.data$minor.copy.number.inTumour[ascat.data$seg_no == first.seg]
    
    for(seg in unique(ascat.data$seg_no[ascat.data$Chromosome == chrom])){
      
      if(ascat.data$total.copy.number.inTumour[ascat.data$seg_no == seg] == nT & 
         ascat.data$minor.copy.number.inTumour[ascat.data$seg_no == seg] == nB){
        
        end <- ascat.data$chromEnd[ascat.data$seg_no == first.seg]
        
      }
      if(ascat.data$total.copy.number.inTumour[ascat.data$seg_no == seg] != nT | 
         ascat.data$minor.copy.number.inTumour[ascat.data$seg_no == seg] != nB){
        
        end <- ascat.data$chromEnd[ascat.data$seg_no == seg]
        
        df_total <- rbind(df_total,data.frame(seg,chrom,start,end,2,1,nT,nB))
        
        start <- ascat.data$chromStart[ascat.data$seg_no == seg]
        nT <- ascat.data$total.copy.number.inTumour[ascat.data$seg_no == seg]
        nB <- ascat.data$minor.copy.number.inTumour[ascat.data$seg_no == seg]
        
      }
    }  
  } 
  names(df_total) <- names(ascat.data)
  return(df_total)
}

summarize_LOH <- function(LOH_seg_file,sample_name,CNVsigs=F){
  print("summarizing LOH")
  ascat.data <- read.table(LOH_seg_file,sep="\t",header=TRUE)
  
  #count number of LOH segments
  LOH_table <- ascatToHRDLOH(ascat.data=ascat.data,
                             SAMPLE.ID=sample_name,
                             return.loc = T)
  
  #save LOH segments and LOH count to ls   
  LOH_ls <- list(nrow(LOH_table),as.matrix(LOH_table[,c("Chromosome","Start","End","totalCN")]))
  names(LOH_ls) <- c("LOHcount","LOHsegments")
  
  if(CNVsigs=T){
    seg <- ascat.data[,c("Chromosome","chromStart","chromEnd","total.copy.number.inTumour","minor.copy.number.inTumour"),] 
    seg$length <- (seg$chromEnd - seg$chromStart)/1000
    seg$Aallele <- seg$total.copy.number.inTumour - seg$minor.copy.number.inTumour
    
    ##Class1
    #"HET" heterozygous segments with copy number of (A > 0, B > 0); 
    seg$BASEclass[seg$minor.copy.number.inTumour > 0 & seg$Aallele > 0] <- "HET"
    #"LOH": segments with LOH with copy number of (A > 0, B = 0);
    seg$BASEclass[seg$minor.copy.number.inTumour == 0 & seg$Aallele > 0] <- "LOH"
    #"HD": segments with homozygous deletions (A = 0, B = 0) 
    seg$BASEclass[seg$minor.copy.number.inTumour == 0 & seg$Aallele == 0] <- "HD"
    
    #Segments were further subclassified into five classes 
    #on the basis of the sum of major and minor alleles (TCN) 
    #TCN = 0 (homozygous deletion); 
    #TCN = 1 (deletion leading to LOH);
    #TCN = 2 (wild type, including copy-neutral LOH); 
    #TCN = 3 or 4 (minor gain); 
    #TCN = 5–8 (moderate gain); and 
    #TCN ≥ 9 (high-level amplification)
    seg$TCNclass[seg$total.copy.number.inTumour == 0] <- "0"
    seg$TCNclass[seg$total.copy.number.inTumour == 1] <- "1"
    seg$TCNclass[seg$total.copy.number.inTumour == 2] <- "2"
    seg$TCNclass[seg$total.copy.number.inTumour == 3 | seg$total.copy.number.inTumour == 4 ] <- "3-4"
    seg$TCNclass[seg$total.copy.number.inTumour >= 5 & seg$total.copy.number.inTumour <= 8] <- "5-8"
    seg$TCNclass[seg$total.copy.number.inTumour >= 9] <- "9+"
    
    # Each of the heterozygous and LOH TCN states were then subclassified into five classes 
    #on basis of the size of their segments: 
    seg$LENGTHclass[seg$length <= 100] <- "0–100kb"
    seg$LENGTHclass[seg$length > 100 & seg$length <= 1000] <- "100kb–1Mb"
    seg$LENGTHclass[seg$length > 1000 & seg$length <= 10000] <- "1Mb–10Mb"
    seg$LENGTHclass[seg$length > 10000 & seg$length <= 40000] <- "10Mb–40Mb"
    seg$LENGTHclass[seg$length > 40000] <- "40Mb+"
    seg$LENGTHclass[seg$length > 1000 & seg$BASEclass == "HD"] <- "1Mb+"
    
    CNV_sigs <- seg %>% group_by(LENGTHclass,TCNclass,BASEclass) %>% tally()
    CNV_sigs$prop <- CNV_sigs$n / sum(CNV_sigs$n)
    
    CNVsigs.class <- fread('~/Documents/GitHub/sigtools_workflow/CNVsigs.class.txt')
    CNVsigs.classifications <- fread('~/Documents/GitHub/sigtools_workflow/COSMIC_v3.3_CN_GRCh37.txt')
    
    CNVsigs.classfied <- left_join(CNVsigs.class,CNV_sigs)
    CNVsigs.classfied[is.na(CNVsigs.classfied)] <- 0
    
    class_prop <- CNVsigs.classfied$prop
    names(class_prop) <- CNVsigs.classfied$Type
  
    #count number of LOH segments
    LOH_table <- ascatToHRDLOH(ascat.data=ascat.data,
                               SAMPLE.ID=sample_name,
                               return.loc = T)
    
    #save LOH segments and LOH count to ls   
    LOH_ls <- list(nrow(LOH_table),as.matrix(LOH_table[,c("Chromosome","Start","End","totalCN")]),class_prop)
    names(LOH_ls) <- c("LOHcount","LOHsegments","CNVclassification")
    
  }


  return(LOH_ls)
}



summarize_SVs <- function(SV_bedpe,tissue){
  
  #start sigtools SV analysis
  SV_bedpe_sigtools <- SV_bedpe
  
  #replace svclass name 
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "DEL"] <- "deletion"
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "DUP"] <- "tandem-duplication" 
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "INV"] <- "inversion" 
  
  #make SV catalog
  SV_catalogue <- bedpeToRearrCatalogue(SV_bedpe_sigtools)
  SV_catalogue_reformat <-  sigTools_formatter(SV_catalogue,sampleName = sample_name)
  
  
  #get organ signature
  SV_sigs <- getOrganSignatures(organ = tissue, 
                                version = "1", 
                                typemut = "rearr", 
                                verbose = FALSE)
  
  #fit organ signature
  SV_fit_sigs <- Fit(catalogues = SV_catalogue_reformat,
                     signatures = SV_sigs,
                     useBootstrap = TRUE,
                     nboot = boots,
                     nparallel = 1)
  
  #extract signature exposures
  SV_exposure_tissue <- as.data.frame(t(SV_fit_sigs$exposures[1,]))
  
  #for CHORD SV analysis
  SV_bedpe_CHORD <- SV_bedpe
  SV_bedpe_CHORD$length <-  SV_bedpe_CHORD$start2 - SV_bedpe_CHORD$start1
  
  
  #make SV list
  SV_ls <- list(nrow(SV_bedpe),SV_catalogue$rearr_catalogue,SV_exposure_tissue,SV_bedpe_CHORD[,c("svclass","length")])

  return(SV_ls)
}

summarize_indels <- function(indel_vcf_location,sample_name,genomeVersion){
  
  
  #sigtools indel classification
  indels_class <- vcfToIndelsClassification(indel_vcf_location,
                                            sampleID =sample_name,
                                            genome.v = genomeVersion
  )
  
  indel_info <- indels_class$count_proportion[c("all.indels","ins","del.mh","del.rep","del.none","del.mh.prop","del.rep.prop","del.none.prop")]
  
  #make indel list
  indel_ls <- list(unlist(indel_info[1]),indel_info,NA)
  return(indel_ls)
}  

summarize_SNVs <- function(snv_vcf_location,genomeVersion,sample_name,tissue){
  #catalog SNVs
  snv_catalogue <- vcfToSNVcatalogue(vcfFilename = snv_vcf_location,
                                     genome.v = genomeVersion)
  
  snv_catalogue_reformat <- sigTools_formatter(input=snv_catalogue,
                                               sampleName=sample_name)
  
  #fit COSMIC v2 signatures, aka 'classic' signatures
  print("fitting classic SNV signatures")
  
  subs_COSMICV2_res <- Fit(catalogues = snv_catalogue_reformat,
                           signatures = COSMIC30_subs_signatures,
                           useBootstrap = TRUE,
                           nboot = boots)
  
  #get organ-specific signatures and fit, a la Degasperi et al 2020
  print("fitting tissue SNV signatures")
  
  SNV_sigs <- getOrganSignatures(organ = tissue, 
                                 version = "1", 
                                 typemut = "subs")
  
  
  SNV_fit_sigs <- Fit(catalogues = snv_catalogue_reformat,
                      signatures = SNV_sigs,
                      useBootstrap = TRUE,
                      nboot = boots)
  
  #add extra step for fitting rare signatures, a la Degasperi et al 2022 
  
  print("fitting rare SNV signatures")
  
  rare_fit <- FitMS(catalogues = snv_catalogue_reformat,
                    organ = tissue,
                    useBootstrap = TRUE,
                    nboot = boots)
  
  #reformat for JSON
  classic_fit <- as.data.frame(t(subs_COSMICV2_res$exposures[1,]))
  
  tissue_fit <- as.data.frame(t(SNV_fit_sigs$exposures[1,]))
  
  rare_fit_df <-  as.data.frame(t(rare_fit$exposures[1,]))
  
  if(length(rare_fit_df) == 1){
    names(rare_fit_df) <- colnames(rare_fit$exposures)
  }
  
  #make list
  SNV_ls <- list(sum(snv_catalogue$catalogue),classic_fit,tissue_fit,rare_fit_df,snv_catalogue_reformat)
  return(SNV_ls)
}  

####arguments####
option_list = list(
  make_option(c("-s", "--sampleName"), type="character", default=NULL, help="sample name", metavar="character"),
  make_option(c("-t", "--tissue"), type="character", default=NULL, help="tissue/organ for HRDetect signatures", metavar="character"),
  make_option(c("-S", "--snvFile"), type="character", default=NULL, help="SNV vcf file", metavar="character"),
  make_option(c("-I", "--indelFile"), type="character", default=NULL, help="indel vcf file", metavar="character"),
  make_option(c("-V", "--SVFile"), type="character", default=NULL, help="Structural variant bed file", metavar="character"),
  make_option(c("-L", "--LOHFile"), type="character", default=NULL, help="LOH file", metavar="character"),
  make_option(c("-b", "--bootstraps"), type="numeric", default=2500, help="number of bootstraps for signature detection", metavar="numeric"),
  make_option(c("-g", "--genomeVersion"), type="character", default="hg38", help="genome version", metavar="character"),
  make_option(c("-i", "--indelCutoff"), type="numeric", default=10, help="minimum number of indels for analysis", metavar="numeric"),
  make_option(c("-c", "--snvCutoff"), type="numeric", default=10, help="minimum number of snvs for analysis", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

boots               <-  opt$bootstraps
genomeVersion       <-  opt$genomeVersion
indelCutoff         <-  opt$indelCutoff
snvCutoff           <-  opt$snvCutoff
sample_name         <-  opt$sampleName
tissue              <-  opt$tissue
snv_vcf_location    <-  opt$snvFile
indel_vcf_location  <-  opt$indelFile
SV_bedpe_location   <-  opt$SVFile
LOH_seg_location    <-  opt$LOHFile

##test
#setwd('/Volumes/')
#boots               <- 2500 
#genomeVersion       <- "hg38"

#indelCutoff         <- 10
#snvCutoff           <- 10
#sample_name         <- "OCT_010434" 
#tissue              <- "Ovary" 
#snv_vcf_location    <- "cgi/scratch/fbeaudry/sigTools_test/TGL62/OCT_010434/OCT_010434_Ov_P_WG.filter.deduped.realigned.recalibrated.mutect2.filtered.VAF.snv.vcf.gz" 
#indel_vcf_location  <- "cgi/scratch/fbeaudry/sigTools_test/TGL62/OCT_010434/OCT_010434_Ov_P_WG.filter.deduped.realigned.recalibrated.mutect2.filtered.VAF.indel.vcf.gz" 
#SV_bedpe_location   <- "cgi/scratch/fbeaudry/sigTools_test/TGL62/OCT_010434/OCT_010434_Ov_P_WG__somatic.somatic_filtered.delly.merged.bedpe"
#LOH_seg_location    <- "cgi/scratch/fbeaudry/sigTools_test/TGL62/OCT_010434/OCT_010434_Ov_P_WG_segments.cna.txt" 

####tissue test####
tissue.catalog <-  c("Biliary", "Bladder", "Bone_SoftTissue", "Breast",  "CNS", "Colorectal", "Esophagus", "Head_neck", "Kidney", "Liver", "Lung", "Lymphoid", "Ovary", "Pancreas", "Prostate", "Skin", "Stomach", "Uterus")

if(tissue %ni% tissue.catalog){
  print("careful: tissue not in tissue-catalog!")
}


####LOH####

LOH_ls <- summarize_LOH(LOH_seg_file=LOH_seg_location,
                        sample_name=sample_name
                        )

####Structural Variants####

SV_bedpe <- try(read.table(SV_bedpe_location,sep = "\t", header = TRUE))

if("try-error" %in% class(SV_bedpe)) {
  
  print("no SVs!")
  SV_ls <- list(0,NA,NA,NA)
  names(SV_ls) <- c("SVcount","SV_sigtools_catalog","SV_sigtools_exposures","SV_CHORD_catalog")
  
} else if(tissue %ni% tissue.catalog){
  
  print("tissue not in tissue-catalog, skipping SV sigs")
  SV_ls <- list(nrow(SV_bedpe),NA,NA,NA)
  names(SV_ls) <- c("SVcount","SV_sigtools_catalog","SV_sigtools_exposures","SV_CHORD_catalog")
  
} else {
  
  print("summarizing SVs")
  SV_ls <- summarize_SVs(SV_bedpe,tissue)
  names(SV_ls) <- c("SVcount","SV_sigtools_catalog","SV_sigtools_exposures","SV_CHORD_catalog")
  
}

####INDEL####

indel_vcf <- try(read.table(indel_vcf_location,comment.char= "#"))

if("try-error" %in% class(indel_vcf) ) {
  
  print("no indels!")
  indel_ls        <- list(0,NA,NA)
  names(indel_ls) <- c("indelCount","indel_sigtools_catalog","indel_CHORD_catalog")
  
} else if(nrow(indel_vcf) < indelCutoff){
  
  print("low indels!")
  indel_ls        <- list(nrow(indel_vcf),NA,NA)
  names(indel_ls) <- c("indelCount","indel_sigtools_catalog","indel_CHORD_catalog")
  
} else {
  
  print("summarizing indels")
  indel_ls        <- summarize_indels(indel_vcf_location,sample_name,genomeVersion)
  names(indel_ls) <- c("indelCount","indel_sigtools_catalog","indel_CHORD_catalog")
  
}

#rename indel file for HRDetect input
names(indel_vcf_location)[1] <- sample_name

####SNV####

snv_vcf <- try(read.table(snv_vcf_location,comment.char= "#"))

if("try-error" %in% class(snv_vcf)  ) {
  
  print("no SNVs!")
  SNV_ls        <- list(0,NA,NA,NA,NA)
  names(SNV_ls) <- c("SNVCount","classic_sigs","tissue_sigs","rare_sigs","SNV_catalog")
  
} else if( nrow(snv_vcf) < snvCutoff){
  
  print("low SNVs!")
  SNV_ls        <- list(nrow(snv_vcf),NA,NA,NA,NA)
  names(SNV_ls) <- c("SNVCount","classic_sigs","tissue_sigs","rare_sigs","SNV_catalog")
  
} else if(tissue %ni% tissue.catalog){
  
  print("tissue not in tissue-catalog, skipping SNV sigs")
  SNV_ls        <- list(nrow(snv_vcf),NA,NA,NA,NA)
  names(SNV_ls) <- c("SNVCount","classic_sigs","tissue_sigs","rare_sigs","SNV_catalog")
  
}  else {
  
  print("summarizing SNVs")
  SNV_ls        <- summarize_SNVs(snv_vcf_location,genomeVersion,sample_name,tissue)
  names(SNV_ls) <- c("SNVCount","classic_sigs","tissue_sigs","rare_sigs","SNV_catalog")
  
}
  
####HRD tests####
HR_calls <- list(NA,NA)
names(HR_calls) <- c("HRDetect","CHORD")

if("try-error" %in% class(snv_vcf) | "try-error" %in% class(indel_vcf) | "try-error" %in% class(SV_bedpe)) {

  print("some data missing, no HRD call!")

  
} else if(nrow(snv_vcf) < snvCutoff | nrow(indel_vcf) < indelCutoff ){

  
  print("some calls too few, no HRD call!")

} else if(tissue %ni% tissue.catalog){
  
  print("tissue not in tissue-catalog, skipping HRD call")
    
}  else {
  print("Performing HRD tests")
  
  #make HRDetect input matrix
  input_matrix <- matrix(NA, nrow = 1, ncol = 6,
                         dimnames = list(
                          sample_name,
                          c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
                        ))
  
  #HRDetect only needs count of LOH segments
  input_matrix[sample_name,"hrd"] <- unlist(LOH_ls[1])
  
  SV_sigtools_catalog <-  sigTools_formatter(SV_ls["SV_sigtools_catalog"],sampleName = sample_name)
  SNV_catalog         <-  sigTools_formatter(SNV_ls["SNV_catalog"],sampleName = sample_name)
  
  #call HRDetect
  HRDetect_res <- HRDetect_pipeline(  SNV_signature_version   = "RefSigv1",
                                      organ                   = tissue,
                                      genome.v                = genomeVersion,
                                      bootstrapHRDetectScores = TRUE,
                                      nbootFit                = boots,
                                      data_matrix             = input_matrix,
                                      SV_catalogues           = SV_sigtools_catalog,
                                      SNV_catalogues          = SNV_catalog,
                                      Indels_vcf_files        = indel_vcf_location
                                    )
  
  #catalog CHORD
  SV_CHORD_catalog <-  as.data.frame(SV_ls["SV_CHORD_catalog"])
  names(SV_CHORD_catalog) <- c("sv_type","sv_len")
  
  CHORD_contexts <- extractSigsChord(
                                      vcf.snv    = snv_vcf_location,
                                      vcf.indel  = indel_vcf_location,
                                      df.sv      = SV_CHORD_catalog,
                                      ref.genome = BSgenome.Hsapiens.UCSC.hg38,
                                      verbose    = F
                                    )     
  
  #call CHORD
  CHORD_call <- chordPredict(CHORD_contexts, do.bootstrap=T, verbose=F,bootstrap.iters = boots)
  
  #assemble list
  HRD_in <- as.data.frame(t(HRDetect_res$data_matrix[1,]))
  
  quantiles <- as.data.frame(t(HRDetect_res$q_5_50_95[1,]))
  colnames(quantiles) <- c("HRD_low_quant","HRD_median","HRD_top_quant")
  
  HRD_point <- HRDetect_res$hrdetect_output[8]
  
  modelWeights <- as.data.frame(t(HRDetect_res$hrdetect_output[1,c(2:7)]))
  
  HRDetect_call <- list(HRD_point,quantiles,modelWeights,HRD_in)
  names(HRDetect_call) <- c("point","quantiles","modelWeights","HRD_in")

  HR_calls <- list(HRDetect_call,CHORD_call)
  names(HR_calls) <- c("HRDetect","CHORD")
  
  #add CHORD catalogs
  SV_ls$SV_CHORD_catalog <- as.data.frame(t(unlist(CHORD_contexts[1,c(127:145)])))
  indel_ls$indel_CHORD_catalog <- as.data.frame(t(unlist(CHORD_contexts[1,c(97:126)])))

}   

#reformat catalogs
if(!is.na(SNV_ls$SNV_catalog)){
  SNV_ls$SNV_catalog <- reformat_toJSON(SNV_ls$SNV_catalog)
}

if(!is.na(SV_ls$SV_sigtools_catalog)){
  SV_ls$SV_sigtools_catalog <- reformat_toJSON(SV_ls$SV_sigtools_catalog)
}

if(is.na(HR_calls$CHORD)){
  SV_ls$SV_CHORD_catalog <- NA
}

#stick all the results together
all.results <- list(sample_name,LOH_ls,SV_ls,indel_ls,SNV_ls,HR_calls)
names(all.results) <- c("Sample","LOH","SV","indel","SNV","HRD")

#conver to JSON and write
ListJSON <- jsonlite::toJSON(all.results,pretty=TRUE,auto_unbox=TRUE)

write(ListJSON,file = paste(sample_name,".signatures.json",sep=""))


