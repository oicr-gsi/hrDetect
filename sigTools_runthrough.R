#! /usr/bin/env Rscript

library(optparse)

# Calling call_hrdetect.R
basedir <- paste(Sys.getenv(c("BASE_DIR")), sep='/')
source(paste0(basedir, "/call_hrdetect.R"))

option_list = list(
  make_option(c("-s", "--sampleName"), type="character", default=NULL, help="sample name", metavar="character"),
  make_option(c("-r", "--SVrefSigs"), type="character", default=NULL, help="structural variant reference signatures", metavar="character"),
  make_option(c("-R", "--SNVrefSigs"), type="character", default=NULL, help="single nucleotide polymorphism reference signatures", metavar="character"),
  make_option(c("-S", "--snvFile"), type="character", default=NULL, help="SNV vcf file", metavar="character"),
  make_option(c("-I", "--indelFile"), type="character", default=NULL, help="indel vcf file", metavar="character"),
  make_option(c("-V", "--SVFile"), type="character", default=NULL, help="Structural variant file", metavar="character"),
  make_option(c("-L", "--LOHFile"), type="character", default=NULL, help="LOH file", metavar="character"),
  make_option(c("-b", "--bootstraps"), type="numeric", default=2500, help="number of bootstraps for signature detection", metavar="numeric"),
  make_option(c("-g", "--genomeVersion"), type="character", default="hg38", help="genome version", metavar="character"),
  make_option(c("-i", "--indelCutoff"), type="numeric", default=10, help="minimum number of indels for analysis", metavar="numeric"),
  make_option(c("-c", "--snvCutoff"), type="numeric", default=10, help="minimum number of snvs for analysis", metavar="numeric"),
  make_option(c("-v", "--svCutoff"), type="numeric", default=10, help="minimum number of SVs for analysis", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

boots               <-  opt$bootstraps
genomeVersion       <-  opt$genomeVersion
indelCutoff         <-  opt$indelCutoff
snvCutoff           <-  opt$snvCutoff
svCutoff            <-  opt$svCutoff
sample_name         <-  opt$sampleName
SVrefSigs           <-  opt$SVrefSigs
SNVrefSigs          <-  opt$SNVrefSigs
snv_vcf_location    <-  opt$snvFile
indel_vcf_location  <-  opt$indelFile
SV_vcf_location     <-  opt$SVFile
LOH_seg_file        <-  opt$LOHFile



## import signatures
  SV_sigs  <- read.table(SVrefSigs, sep = "\t", header = TRUE)
  SNV_sigs <- read.table(SNVrefSigs, sep = "\t", header = TRUE)

#### PREPROCESS ALL FOUR DATA TYPES ####
  
missing_data = FALSE
  
#### pre-process SEGMENT data ####
  cat('pre-processing LOH data\n')
  seg.data <- read.table(LOH_seg_file, sep="\t", header=TRUE)

#### pre-process In/Del data ####
  cat('pre-processing in/del data\n')
  indel_vcf <- try(read.table(indel_vcf_location, sep = "\t", header = TRUE))
  
  if("try-error" %in% class(indel_vcf)){
    ID.catalog.JSON <- jsonlite::toJSON(list("QC"="FAIL","count"=0,"Results"="NA"),pretty=TRUE,auto_unbox=TRUE)
    missing_data = TRUE
  }else{
    indels_class    <- vcfToIndelsClassification(indel_vcf_location, sample_name, genome.v = genomeVersion)
    ID.catalog.JSON <- jsonlite::toJSON(list("QC"="PASS","count"=nrow(indel_vcf),"Results"=indels_class$count_proportion),pretty=TRUE,auto_unbox=TRUE)
    if( nrow(indel_vcf) < indelCutoff ){missing_data = TRUE}
    
  }

  write(ID.catalog.JSON, file = paste0(sample_name,".catalog.ID.json"))

#### pre-process Structural Variant data ####
  cat('pre-processing SV data\n')
  
  SV_vcf        <- try(read.table(SV_vcf_location, sep = "\t", header = TRUE))
  
  if("try-error" %in% class(SV_vcf)){
    SV.JSON <- jsonlite::toJSON(list("QC"="FAIL","count"=0,"catalog"="NA","exposures"="NA"),pretty=TRUE,auto_unbox=TRUE)  
    missing_data = TRUE
    
  }else{  
    SV_catalogue  <- SV_vcf_cataloger(SV_vcf_location, sample_name)
    SV_ls     <- summarize_SVs(SV_catalogue, SV_sigs, sample_name, nboot = 2000)
    
    SV_ls.preJSON <- list("QC"="PASS",
                          "count"= nrow(SV_vcf),
                          "catalog"=as.data.frame(t(SV_ls$SV_sigtools_catalog[1])),
                          "exposures"=SV_ls$SV_sigtools_exposures)
    SV.JSON <- jsonlite::toJSON(SV_ls.preJSON,pretty=TRUE,auto_unbox=TRUE)
    if( nrow(SV_vcf) < svCutoff ){missing_data = TRUE}
    
  }
  write(SV.JSON, file = paste0(sample_name,".exposures.SV.json"))

#### pre-process SNV data ####
  cat('pre-processing SNV data\n')
  
  snv_df    <- try(read.table(snv_vcf_location, sep = "\t", header = FALSE))
  
  if("try-error" %in% class(snv_df)){
    
    SBS.JSON <- jsonlite::toJSON(list("QC"="FAIL","count"=0,"catalog"="NA","exposures"="NA"),pretty=TRUE,auto_unbox=TRUE)  
    missing_data = TRUE
    
  }else{  
    names(snv_df)[c(1,2,4,5)] <- c("Chromosome", "Start_Position", "Reference_Allele", "Allele")
    snv_df$Tumor_Sample_Barcode <- sample_name
    
    snv_catalogue <- mut.to.sigs.input(mut.ref=snv_df,
                                       sample.id="Tumor_Sample_Barcode",
                                       chr="Chromosome",
                                       pos="Start_Position",
                                       ref="Reference_Allele",
                                       alt="Allele",
                                       bsg=BSgenome.Hsapiens.UCSC.hg38)
    
    SNV_ls    <- summarize_SNVs_dcSigs(snv_df, SNV_sigs, sample_name, subsample_proportion=1)
    SNV_ls.preJSON <- list("QC"="PASS",
                           "count"= nrow(snv_df),
                          "catalog"=snv_catalogue,
                          "exposures"=SNV_ls
                          )
    
    SBS.JSON <- jsonlite::toJSON(SNV_ls.preJSON,pretty=TRUE,auto_unbox=TRUE)
    if( nrow(snv_df) < snvCutoff ){missing_data = TRUE}
  }
  write(SBS.JSON, file = paste0(sample_name,".exposures.SBS.json"))

#### run HRDetect ####
  
if(missing_data == TRUE) {
  
  cat("some data missing, no HRD call!\n")
  
  hrdetect.json <- jsonlite::toJSON(
                  list("sample_name"=sample_name,
                        "QC"="FAIL",
                        "hrdetect_call"="NA"
  ),pretty=TRUE,auto_unbox=TRUE)
  
} else {
    
  quantiles <-
    call_hrdetect(boots       ,
                sample_name   ,
                SV_sigs       ,
                SNV_sigs      ,
                snv_df        ,
                indels_class  ,
                SV_catalogue  ,
                seg.data        
   )
  hrdetect_call <- list("sample_name"=sample_name,
                        "QC"="PASS",
                        "hrdetect_call"=quantiles
  )

  hrdetect.json <- jsonlite::toJSON(hrdetect_call,pretty=TRUE,auto_unbox=TRUE)
  
}

write(hrdetect.json, file = paste0(sample_name,".signatures.json"))
  