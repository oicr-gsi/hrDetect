##version 1.2

####packages####
#install_github('Nik-Zainal-Group/signature.tools.lib',ref='v2.1.2')
library(signature.tools.lib)
library(optparse)

####functions####

sigTools_formatter <- function(input,sampleName){
  cat_list <- list()
  names(input)[1] <- "catalogue"
  colnames(input$catalogue) <- sampleName
  cat_list[[1]] <- input$catalogue
  catalogues <- do.call(cbind,cat_list)
  return(catalogues)
}

setColNames <- function (object = nm, nm) {
  colnames(object) <- nm
  object
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
  make_option(c("-i", "--indelCutoff"), type="numeric", default=10, help="minimum number of indels for analysis", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

sample_name <- opt$sampleName
tissue <-  opt$tissue
snvFile_loc  <-  opt$snvFile
indel_vcf_file <-  opt$indelFile
SV_bedpe_file <-  opt$SVFile
LOH_seg_file <-  opt$LOHFile
boots <-  opt$bootstraps
genomeVersion <-  opt$genomeVersion
indelCutoff <-  opt$indelCutoff

##test
#sample_name <-  "PANX_1309" 
#tissue <- "Pancreas" 
#snvFile_loc  <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309_Lv_M_WG_100-PM-033_LCM4.filter.deduped.realigned.recalibrated.mutect2.filtered.VAF.SNP.vcf.gz" 
#indel_vcf_file <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309_Lv_M_WG_100-PM-033_LCM4.filter.deduped.realigned.recalibrated.mutect2.filtered.VAF.indel.vcf.gz" 
#SV_bedpe_file <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309_Lv_M_WG_100-PM-033_LCM4_somatic.somatic_filtered.delly.merged.bedpe"
#LOH_seg_file <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309_Lv_M_WG_100-PM-033_LCM4_segments.cna.txt" 
#boots <- 2500 
#genomeVersion <- "hg38"
#indelCutoff <- 10

####Import files####

#HRDetect_pipeline() will throw error if there are no indels
t <- try(read.table(indel_vcf_file,comment.char= "#"))
if("try-error" %in% class(t) | length(t) < indelCutoff) {
  
  for( fileType in c("hrd","model","sigs")){
    write.table(
      c("no indels after filtering"),
      file = paste(sample_name,".sigtools.",fileType,".txt",sep=""),append = F, quote = FALSE, row.names = FALSE, col.names = FALSE
    )
  }
  
}else{
  
  ####Take in files####
  
  ##SNV##
  
  snv_catalogue <- vcfToSNVcatalogue(vcfFilename = snvFile_loc,genome.v = genomeVersion)
  
  snv_catalogue_reformat <- sigTools_formatter(input=snv_catalogue,sampleName=sample_name)
  
  ##structural variants##
  
  SV_bedpe <- read.table(SV_bedpe_file,
                         sep = "\t",header = TRUE,
                         stringsAsFactors = FALSE,check.names = FALSE)
  
  #replace svclass name 
  SV_bedpe$svclass[SV_bedpe$svclass == "DEL"] <- "deletion"
  SV_bedpe$svclass[SV_bedpe$svclass == "DUP"] <- "tandem-duplication" 
  SV_bedpe$svclass[SV_bedpe$svclass == "INV"] <- "inversion" 
  
  SV_catalogue <- bedpeToRearrCatalogue(SV_bedpe)
  
  SV_catalogue_reformat <- sigTools_formatter(input=SV_catalogue,sampleName=sample_name)
  
  ##LOH##
  
  ascat.data <- read.table(LOH_seg_file,sep="\t",header=TRUE)
  
  hrd_index <- ascatToHRDLOH(ascat.data=ascat.data,SAMPLE.ID=sample_name)
  
  ##INDEL##
  
  names(indel_vcf_file)[1] <- sample_name
  
  ####HRD test####
  
  col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  
  input_matrix <- matrix(NA,nrow = 1,
                         ncol = length(col_hrdetect),
                         dimnames = list(sample_name,col_hrdetect))
  
  input_matrix[sample_name,"hrd"] <- hrd_index
  
  HRDetect_res <- HRDetect_pipeline(data_matrix=input_matrix,
                                      bootstrapHRDetectScores=TRUE,
                                      SV_catalogues=SV_catalogue_reformat,
                                      SNV_catalogues=snv_catalogue_reformat,
                                      organ=tissue,
                                      Indels_vcf_files=indel_vcf_file,
                                      genome.v = genomeVersion,
                                      nbootFit=boots
                                      )
    
  ####output####
    
  quantiles <- HRDetect_res$q_5_50_95
  names(quantiles) <- c("HRD_low_quant","HRD_median","HRD_top_quant")
    
  write.table(
    c(
        "HRD_point"=HRDetect_res$hrdetect_output[8],
        quantiles
      ),
      file = paste(sample_name,".sigtools.hrd.txt",sep=""),
      append = F, quote = FALSE, sep = "\t", 
      eol = "\n", na = "NA",dec = ".", row.names = TRUE, 
      col.names = FALSE
    )
    
  write.table(
      rbind(
      "value"=HRDetect_res$data_matrix,
      "BLUP"=HRDetect_res$hrdetect_output[c(2:7)]
      ),
      file = paste(sample_name,".sigtools.model.txt",sep=""),
      append = F, quote = FALSE, sep = "\t", 
      eol = "\n", na = "NA",dec = ".", row.names = TRUE, 
      col.names = TRUE
  )
  
  ## signature objects will be empty if too few SNPs or SVs
  if(length(HRDetect_res$exposures_rearr)>0){
    
    if(length(HRDetect_res$exposures_subs)>0){
      
      write.table(
        rbind.data.frame(
          cbind.data.frame(
            setColNames(HRDetect_res$exposures_subs, "sig_weight_norm"),
            setColNames(HRDetect_res$exposures_subs/sum(HRDetect_res$exposures_subs), "sig_weight_rel"),
            setColNames(HRDetect_res$exposures_subs/sum(HRDetect_res$SNV_catalogues), "sig_weight_rel_adj"),
            "sig_type"="SNV"
          ),
          cbind.data.frame(
            setColNames(HRDetect_res$exposures_rearr, "sig_weight_norm"),
            setColNames(HRDetect_res$exposures_rearr/sum(HRDetect_res$exposures_rearr), "sig_weight_rel"),
            setColNames(HRDetect_res$exposures_rearr/sum(HRDetect_res$SV_catalogues), "sig_weight_rel_adj"),
            "sig_type"="structural"
          )
        ),
        file = paste(sample_name,".sigtools.sigs.txt",sep=""),
        append = F, quote = FALSE, sep = "\t", 
        eol = "\n", na = "NA",dec = ".", row.names = TRUE, 
        col.names = TRUE
      )
      
    }
    
    if(length(HRDetect_res$exposures_subs)==0){
      
      write.table(
          cbind.data.frame(
            setColNames(HRDetect_res$exposures_rearr, "sig_weight_norm"),
            setColNames(HRDetect_res$exposures_rearr/sum(HRDetect_res$exposures_rearr), "sig_weight_rel"),
            setColNames(HRDetect_res$exposures_rearr/sum(HRDetect_res$SV_catalogues), "sig_weight_rel_adj"),
            "sig_type"="structural"
          ),
        file = paste(sample_name,".sigtools.sigs.txt",sep=""),
        append = F, quote = FALSE, sep = "\t", 
        eol = "\n", na = "NA",dec = ".", row.names = TRUE, 
        col.names = TRUE
      )
      
    }  
  }
    
  if(length(HRDetect_res$exposures_rearr)==0){  
    
      write.table(
          cbind.data.frame(
            setColNames(HRDetect_res$exposures_subs, "sig_weight_norm"),
            setColNames(HRDetect_res$exposures_subs/sum(HRDetect_res$exposures_subs), "sig_weight_rel"),
            setColNames(HRDetect_res$exposures_subs/sum(HRDetect_res$SNV_catalogues), "sig_weight_rel_adj"),
            "sig_type"="SNV"
          ),
        file = paste(sample_name,".sigtools.sigs.txt",sep=""),
        append = F, quote = FALSE, sep = "\t", 
        eol = "\n", na = "NA",dec = ".", row.names = TRUE, 
        col.names = TRUE
      )
    
  }
    
} #end of if(few indels)

