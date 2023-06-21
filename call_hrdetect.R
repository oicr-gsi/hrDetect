#! /usr/bin/env Rscript
#devtools::install_github('Nik-Zainal-Group/signature.tools.lib',ref='v2.1.2')

library(signature.tools.lib)
library(BSgenome.Hsapiens.UCSC.hg38)
library(jsonlite)

'%ni%' <- function(x,y)!('%in%'(x,y))

call_hrdetect <- function(boots               ,
                          genomeVersion       ,
                          indelCutoff         ,
                          snvCutoff           ,
                          sample_name         ,
                          SVrefSigs           ,
                          snv_vcf_location    ,
                          indel_vcf_location  ,
                          SV_bedpe_location   ,
                          LOH_seg_location    
){
  
  ####LOH####
  
  LOH_ls <- summarize_LOH(LOH_seg_file=LOH_seg_location,
                          sample_name=sample_name
  )
  
  ####Structural Variants####
  
  SV_bedpe <- try(read.table(SV_bedpe_location,sep = "\t", header = TRUE))
  
  if("try-error" %in% class(SV_bedpe)) {
    
    print("no SVs!")
    SV_ls <- list(0,NA,NA)
    names(SV_ls) <- c("SVcount","SV_sigtools_catalog","SV_sigtools_exposures")
    
  } else {
    
    print("summarizing SVs")
    SV_ls <- summarize_SVs(SV_bedpe,SVrefSigs)
    names(SV_ls) <- c("SVcount","SV_sigtools_catalog","SV_sigtools_exposures")
    
  }
  
  ####INDEL####
  
  indel_vcf <- try(read.table(indel_vcf_location,comment.char= "#"))
  
  if("try-error" %in% class(indel_vcf) ) {
    
    print("no indels!")
    indel_ls        <- list(0,NA)
    names(indel_ls) <- c("indelCount","indel_sigtools_catalog")
    
  } else if(nrow(indel_vcf) < indelCutoff){
    
    print("low indels!")
    indel_ls        <- list(nrow(indel_vcf),NA)
    names(indel_ls) <- c("indelCount","indel_sigtools_catalog")
    
  } else {
    
    print("summarizing indels")
    indel_ls        <- summarize_indels(indel_vcf_location,sample_name,genomeVersion)
    names(indel_ls) <- c("indelCount","indel_sigtools_catalog")
    
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
    
  }  else {
    
    print("summarizing SNVs")
    SNV_ls        <- summarize_SNVs(snv_vcf_location=snv_vcf_location,
                                    genomeVersion=genomeVersion,
                                    sample_name=sample_name
    )
    names(SNV_ls) <- c("SNVCount","classic_sigs","SNV_catalog")
    
  }
  
  ####HRD tests####
  HR_calls <- list(NA)
  names(HR_calls) <- c("HRDetect")
  
  if("try-error" %in% class(snv_vcf) | "try-error" %in% class(indel_vcf) | "try-error" %in% class(SV_bedpe)) {
    
    print("some data missing, no HRD call!")
    
  } else if(nrow(snv_vcf) < snvCutoff | nrow(indel_vcf) < indelCutoff ){
    
    print("some calls too few, no HRD call!")
    
  } else {
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
                                        genome.v                = genomeVersion,
                                        bootstrapHRDetectScores = TRUE,
                                        nbootFit                = boots,
                                        data_matrix             = input_matrix,
                                        SV_catalogues           = SV_sigtools_catalog,
                                        SNV_catalogues          = SNV_catalog,
                                        Indels_vcf_files        = indel_vcf_location
    )
    
    #assemble list
    HRD_in <- as.data.frame(t(HRDetect_res$data_matrix[1,]))
    
    quantiles <- as.data.frame(t(HRDetect_res$q_5_50_95[1,]))
    colnames(quantiles) <- c("HRD_low_quant","HRD_median","HRD_top_quant")
    
    HRD_point <- HRDetect_res$hrdetect_output[8]
    
    modelWeights <- as.data.frame(t(HRDetect_res$hrdetect_output[1,c(2:7)]))
    
    HRDetect_call <- list(HRD_point,quantiles,modelWeights,HRD_in)
    names(HRDetect_call) <- c("point","quantiles","modelWeights","HRD_in")
    
    
  }   
  
  #reformat catalogs
  if(!is.na(SNV_ls$SNV_catalog)){
    SNV_ls$SNV_catalog <- reformat_toJSON(SNV_ls$SNV_catalog)
  }
  
  if(!is.na(SV_ls$SV_sigtools_catalog)){
    SV_ls$SV_sigtools_catalog <- reformat_toJSON(SV_ls$SV_sigtools_catalog)
  }
  
  #stick all the results together
  all.results <- list(sample_name,LOH_ls,SV_ls,indel_ls,SNV_ls,HR_calls)
  names(all.results) <- c("Sample","LOH","SV","indel","SNV","HRD")
  
  #conver to JSON and write
  ListJSON <- jsonlite::toJSON(all.results,pretty=TRUE,auto_unbox=TRUE)
  
  write(ListJSON,file = paste(sample_name,".signatures.json",sep=""))
}

reformat_toJSON <- function(inCol){
  catalog <-   as.data.frame(t(inCol[,1])) 
  names(catalog) <- rownames(inCol) 
  return(catalog)
}

segtoASCAT <- function(seg.data, input_format){
  col_names <- c( "seg_no", "Chromosome", "chromStart", "chromEnd", 
                  "total.copy.number.inNormal","minor.copy.number.inNormal",
                  "total.copy.number.inTumour","minor.copy.number.inTumour" )
 
   if(input_format == "SEQUENZA"){
    # tail -n +2 ~{segFile} | \
    # awk 'split($1,a,"\"") split(a[2],b,"chr") \
    # c(b[2],$2,$3,2,1,$10,$12)            {print NR"\t"}' 
  }else if(input_format == "PURPLE"){
    
  }
  
}

setColNames <- function (object = nm, nm) {
  
  colnames(object) <- nm
  object
  
}

sigTools_formatter <- function(input,sampleName){
  
  cat_list <- list()
  names(input)[1] <- "catalogue"
  colnames(input$catalogue) <- sampleName
  cat_list[[1]] <- input$catalogue
  catalogues <- do.call(cbind,cat_list)
  return(catalogues)
  
}

summarize_indels <- function(indel_vcf_location, sample_name, genomeVersion){
  
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

summarize_LOH <- function(LOH_seg_file, sample_name){
  
  print("summarizing LOH")
  seg.data <- read.table(LOH_seg_file,sep="\t",header=TRUE)
  
  LOH.data <- seg.data[seg.data$minorAlleleCopyNumber < 0.5 & seg.data$majorAlleleCopyNumber >= 0.5]
  
  LOH.data.big <- LOH.data[LOH.data$end - LOH.data$start > 1.5e+07,]
  # I don't love this length based filter, try LOH grouped by arm-level
  # and then filtered by length - the issue is currently the length of segment
  # will impact filtering but it should not
  
  #save LOH count to ls   
  LOH_ls <- list(nrow(LOH.data.big))
  names(LOH_ls) <- c("LOHcount")
  
  return(LOH_ls)
  
}

summarize_SNVs <- function(snv_vcf_location, genomeVersion, sample_name){
  
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
                           giniThresholdScaling = 2, 
                           threshold_percent = 1,
                           nboot = boots)
  
  #reformat for JSON
  classic_fit <- as.data.frame(t(subs_COSMICV2_res$exposures[1,]))
  
  #make list
  SNV_ls <- list(sum(snv_catalogue$catalogue),classic_fit,snv_catalogue_reformat)
  return(SNV_ls)
  
} 

summarize_SVs <- function(SV_bedpe,SVrefSigs){
  
  #start sigtools SV analysis
  SV_bedpe_sigtools <- SV_bedpe
  
  #replace svclass name 
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "DEL"] <- "deletion"
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "DUP"] <- "tandem-duplication" 
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "INV"] <- "inversion" 
  
  #make SV catalog
  SV_catalogue <- bedpeToRearrCatalogue(SV_bedpe_sigtools)
  SV_catalogue_reformat <-  sigTools_formatter(SV_catalogue,sampleName = sample_name)
  
  SV_sigs <- read.table(SVrefSigs,sep = "\t", header = TRUE)
  
  #fit organ signature
  SV_fit_sigs_all <- Fit(catalogues = SV_catalogue_reformat,
                         signatures = SV_sigs,
                         useBootstrap = TRUE,
                         nboot = boots,
                         giniThresholdScaling=2 ,
                         threshold_percent=1,
                         nparallel = 1)
  
  #extract signature exposures
  SV_exposure <- as.data.frame(t(SV_fit_sigs_all$exposures[1,]))
  
  
  #make SV list
  SV_ls <- list(nrow(SV_bedpe),SV_catalogue$rearr_catalogue,SV_exposure)
  
  return(SV_ls)
}


