#! /usr/bin/env Rscript
#devtools::install_github('Nik-Zainal-Group/signature.tools.lib',ref='v2.1.2')
#options(timeout=1000)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library(signature.tools.lib)
library(BSgenome.Hsapiens.UCSC.hg38)
library(jsonlite)
library(deconstructSigs)
library(dplyr)

'%ni%' <- function(x,y)!('%in%'(x,y))

call_hrdetect <- function(boots               ,
                          genomeVersion       ,
                          sample_name         ,
                          SVrefSigs           ,
                          SNVrefSigs          ,
                          snv_vcf_location    ,
                          indel_vcf_location  ,
                          SV_bedpe_location   ,
                          LOH_seg_file    
){
  print('pre-processing data')
  
  # pre-process data
  seg.data      <- read.table(LOH_seg_file,sep="\t",header=TRUE)
  indels_class  <- vcfToIndelsClassification(indel_vcf_location, sample_name, genome.v = genomeVersion)
  snv_catalogue <- vcfToSNVcatalogue(vcfFilename = snv_vcf_location, genome.v = genomeVersion)
  
  df_snv <- try(read.table(snv_vcf_location,sep = "\t", header = FALSE))
  names(df_snv)[c(1,2,4,5)] <- c("Chromosome","Start_Position","Reference_Allele","Allele")

  #replace svclass name 
  SV_bedpe$svclass[SV_bedpe$svclass == "DEL"] <- "deletion"
  SV_bedpe$svclass[SV_bedpe$svclass == "DUP"] <- "tandem-duplication" 
  SV_bedpe$svclass[SV_bedpe$svclass == "INV"] <- "inversion" 
  SV_catalogue <- bedpeToRearrCatalogue(SV_bedpe)
  
  SV_sigs  <- read.table(SVrefSigs, sep = "\t", header = TRUE)
  SNV_sigs <- read.table(SNVrefSigs, sep = "\t", header = TRUE, row.names=1)
   
  hrdetect_boots <- c("del.mh.prop"  =numeric(),"SNV3"  =numeric(),"SV3"  =numeric(),"SV5"  =numeric(),"hrd"  =numeric(),"SNV8"  =numeric(),"intercept"  =numeric(),
                      "del.mh.prop.w"=numeric(),"SNV3.w"=numeric(),"SV3.w"=numeric(),"SV5.w"=numeric(),"hrd.w"=numeric(),"SNV8.w"=numeric(),"Probability"=numeric())
  
  cat('starting bootstrapping')
  for(boot in c(1:boots)){
    if(boot %% 10 == 0){ cat(boot)}else{ cat('.')  }
    
    ## Summarize data
    LOH_ls    <- summarize_LOH( seg.data)
    indel_ls  <- summarize_indels( indels_class )
    SV_ls     <- summarize_SVs(SV_catalogue, SV_sigs, sample_name)
    #SNV_ls   <- summarize_SNVs(snv_catalogue, SNV_sigs, sample_name)
    SNV_ls    <- dcSigs(df_snv, sample_name)
    
    ## Perform test
    input_matrix <- matrix(NA, nrow = 1, ncol = 6,
                            dimnames = list(
                             sample_name,
                             c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
                            )
                           )
    
    input_matrix[sample_name,"del.mh.prop"] <- indel_ls$sim_proportions["del.mh.prop"]
    #input_matrix[sample_name,"SNV3"] <- SNV_ls$classic_sigs$SBS3
    input_matrix[sample_name,"SNV3"] <- SNV_ls$Signature.3
    
    input_matrix[sample_name,"SV3"]  <- SV_ls$SV_sigtools_exposures$RefSigR3
    input_matrix[sample_name,"SV5"]  <- SV_ls$SV_sigtools_exposures$RefSigR5
    input_matrix[sample_name,"hrd"]  <- unlist(LOH_ls[1])
    # input_matrix[sample_name,"SNV8"] <- SNV_ls$classic_sigs$SBS8 
     input_matrix[sample_name,"SNV8"] <- SNV_ls$Signature.8
    
    suppressMessages(
      HRDetect_res <- HRDetect_pipeline(data_matrix = input_matrix, bootstrapHRDetectScores = FALSE)
    )

    boot_results <- c(HRDetect_res$data_matrix, HRDetect_res$hrdetect_output)
    hrdetect_boots <- rbind.data.frame(hrdetect_boots,boot_results )
    names(hrdetect_boots) <- c(colnames(HRDetect_res$data_matrix), paste0(colnames(HRDetect_res$hrdetect_output),".w"))
    
  }   
  cat('done bootstrapping')
  
  quantiles <- lapply(hrdetect_boots, function (x) quantile(x, probs = c(5, 50, 95)/100))
  
  quantiles <- c( "sample_name"=sample_name, quantiles  )
  quantiles.json <- jsonlite::toJSON(quantiles,pretty=TRUE,auto_unbox=TRUE)
  write(quantiles.json, file = paste0(sample_name,".signatures.json"))
  
 # SBS.signatures <- c("sample_name"=sample_name, SNV_ls$classic_sigs)
  SV.signatures <- c("sample_name"=sample_name, SV_ls$SV_sigtools_exposures)
  
 # SBS.signatures <- jsonlite::toJSON(SBS.signatures,pretty=TRUE,auto_unbox=TRUE)
  SV.signatures <- jsonlite::toJSON(SV.signatures,pretty=TRUE,auto_unbox=TRUE)
  
#  write(SBS.signatures, file = paste0(sample_name,".exposures.SBS.json"))
  write(SV.signatures, file = paste0(sample_name,".exposures.SV.json"))
  
  write.table(SV_ls$SV_sigtools_catalog, file = paste0(sample_name,".catalog.SV.json"), row.names = TRUE)
  write.table(snv_catalogue$catalogue, file = paste0(sample_name,".catalog.SBS.json"), row.names = TRUE)
  write.table(indels_class$count_proportion, file = paste0(sample_name,".catalog.ID.json"), row.names = TRUE)
  
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

sigTools_formatter <- function(input,sampleName){
  
  cat_list <- list()
  names(input)[1] <- "catalogue"
  colnames(input$catalogue) <- sampleName
  cat_list[[1]] <- input$catalogue
  catalogues <- do.call(cbind,cat_list)
  return(catalogues)
  
}

summarize_indels <- function(indels_class){
  
  #sigtools indel classification
  
  indel_info <- indels_class$count_proportion[c("all.indels","ins","del.mh","del.rep","del.none","del.mh.prop","del.rep.prop","del.none.prop")]
  
  
  sim_proportions <- table(sample(1:3 , indel_info$all.indels, replace=TRUE, prob=indel_info[c(6:8)])) / indel_info$all.indels
  names(sim_proportions) <-  c("del.mh.prop","del.rep.prop","del.none.prop")
  
  
  #make indel list
  indel_ls <- list(unlist(indel_info[1]), indel_info, sim_proportions)
  names(indel_ls) <- c("indelCount","indel_sigtools_catalog", "sim_proportions")
  

  return(indel_ls)
}  

summarize_LOH <- function(seg.data){
  library(dplyr)
  
  seg.data$length <- seg.data$end - seg.data$start
  sum_of_segment_lengths <- sum(seg.data$length)
  
  LOH.data <- seg.data[seg.data$minorAlleleCopyNumber < 0.5 & seg.data$majorAlleleCopyNumber >= 0.5,]
  sum_of_LOH_lengths <- sum(LOH.data$length)
  
  LOH_probabilities <- c(sum_of_LOH_lengths/sum_of_segment_lengths, 1-sum_of_LOH_lengths/sum_of_segment_lengths)

  LOH_ls <- table(sample(1:2 , 60, replace=TRUE, prob=LOH_probabilities))[1] 
  # 60 is arbitrary
  
  names(LOH_ls) <- c("LOHcount")
  
  return(LOH_ls)
  
}

summarize_SNVs <- function(snv_catalogue, SNV_sigs, sample_name, giniThresholdScaling=2, threshold_percent=1){
  

  col_names <- c("SNVCount","classic_sigs","SNV_catalog")
  
  snv_catalogue_reformat <- sigTools_formatter(input=snv_catalogue, sampleName=sample_name)
  
  subs_res <- Fit(catalogues = snv_catalogue_reformat,
                           signatures = SNV_sigs,
                           useBootstrap = TRUE, 
                           nboot = 200,
                          exposureFilterType = "giniScaledThreshold",
                          giniThresholdScaling = giniThresholdScaling, 
                          threshold_percent = threshold_percent
                  )
  
  #reformat for JSON
  classic_fit <- as.data.frame(t(subs_res$exposures[1,]))
  
  #make list
  SNV_ls <- list(sum(snv_catalogue$catalogue),classic_fit,snv_catalogue_reformat)
  names(SNV_ls) <- col_names
  
  return(SNV_ls)
  
} 

summarize_SVs <- function(SV_catalogue, SV_sigs, sampleName , giniThresholdScaling=2, threshold_percent=1){
    

  SV_catalogue_reformat <-  sigTools_formatter(SV_catalogue, sampleName = sample_name)
  
  #fit SV signature
  SV_fit_sigs_all <- Fit(catalogues = SV_catalogue_reformat,
                         signatures = SV_sigs,
                         useBootstrap = TRUE,
                         nboot = 200,
                         exposureFilterType = "giniScaledThreshold",
                         giniThresholdScaling = giniThresholdScaling, 
                         threshold_percent = threshold_percent
                        )
  
  #extract signature exposures
  SV_exposure <- as.data.frame(t(SV_fit_sigs_all$exposures[1,]))
  
  #make SV list
  SV_ls <- list(nrow(SV_bedpe),SV_catalogue$rearr_catalogue,SV_exposure)
  
  names(SV_ls) <- c("SVcount","SV_sigtools_catalog","SV_sigtools_exposures")
    
  return(SV_ls)
  
}

dcSigs <- function(df_snv, sample_name){
  
  n_snv <- nrow(df_snv)
  
  # fix formatting
  df_snv$Tumor_Sample_Barcode <- sample_name
  
  df_snv <- dplyr::sample_n(df_snv, round(n_snv*0.8), replace = F)
  
  # read input data
  sigs_input <- mut.to.sigs.input(mut.ref=df_snv,
                                  sample.id="Tumor_Sample_Barcode",
                                  chr="Chromosome",
                                  pos="Start_Position",
                                  ref="Reference_Allele",
                                  alt="Allele",
                                  bsg=BSgenome.Hsapiens.UCSC.hg38)
  
  # calculate a sample
  sample_sigs <- whichSignatures(tumor.ref=sigs_input,
                                 signatures.ref=signatures.nature2013,
                                 sample.id=sample_name,
                                 contexts.needed=TRUE,
                                 tri.counts.method='exome')
  
  
  # get signature weights
  weights <- sample_sigs$weights
  df_nums <- weights * n_snv
  
  return(df_nums)
  
}

simpleEventType <- function(gr) {
  #https://github.com/PapenfussLab/gridss/blob/master/example/simple-event-annotation.R
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
                ifelse(strand(gr) == strand(pgr), "INV",
                       ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
                              ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
                                     "DUP")))))
}

