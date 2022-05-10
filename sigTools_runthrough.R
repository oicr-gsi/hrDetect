##version 1.1

####packages####

library(signature.tools.lib)

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

args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
tissue <- args[2]
snvFile_loc  <- args[3]
indel_vcf_file <- args[4]
SV_bedpe_file <- args[5]
LOH_seg_file <- args[6]
boots <- as.numeric(args[7])
genomeVersion <- args[8]

####Import files####

#HRDetect_pipeline() will throw error if there are no indels, which happens
t <- try(read.table(indel_vcf_file,comment.char= "#"))
if("try-error" %in% class(t)) {
  
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
                                      signature_type=tissue,
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
    
}

