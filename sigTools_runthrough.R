##version 1.3

####packages####
#install_github('Nik-Zainal-Group/signature.tools.lib',ref='v2.1.2')
library(signature.tools.lib)
library(optparse)

library(CHORD)
library(BSgenome.Hsapiens.UCSC.hg38)
library(jsonlite)

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

extend_segments <-  function(ascat.data){
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

sample_name    <-  opt$sampleName
tissue         <-  opt$tissue
snvFile_loc    <-  opt$snvFile
indel_vcf_file <-  opt$indelFile
SV_bedpe_file  <-  opt$SVFile
LOH_seg_file   <-  opt$LOHFile
boots          <-  opt$bootstraps
genomeVersion  <-  opt$genomeVersion
indelCutoff    <-  opt$indelCutoff

##test
#setwd('/Volumes/')
#sample_name <-  "PANX_1309" 
#tissue <- "Pancreas" 
#snvFile_loc  <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309_Lv_M_WG_100-PM-033_LCM4.filter.deduped.realigned.recalibrated.mutect2.filtered.VAF.SNP.vcf.gz" 
#indel_vcf_file <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309_Lv_M_WG_100-PM-033_LCM4.filter.deduped.realigned.recalibrated.mutect2.filtered.VAF.indel.vcf.gz" 
#SV_bedpe_file <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309_Lv_M_WG_100-PM-033_LCM4_somatic.somatic_filtered.delly.merged.bedpe"
#LOH_seg_file <- "cgi/scratch/fbeaudry/sigTools_test/PASS01/PANX_1309/PANX_1309_Lv_M_WG_100-PM-033_LCM4_segments.cna.txt" 
#boots <- 2500 
#genomeVersion <- "hg38"
#indelCutoff <- 10

####LOH####
print("summarizing LOH")
ascat.data <- read.table(LOH_seg_file,sep="\t",header=TRUE)
  
#extend disrupted segments by assuming homozygosity within missing chunk 
#is same on both sides of missing chunk; increases # of LOH segs
ascat.data.ext <- extend_segments(ascat.data)

#count number of LOH segments
LOH_table <- ascatToHRDLOH(ascat.data=ascat.data.ext,
                           SAMPLE.ID=sample_name,
                           return.loc = T)
  
#save LOH segments and LOH count to ls   
LOH_ls <- list(nrow(LOH_table),as.matrix(LOH_table[,c("Chromosome","Start","End","totalCN")]))
names(LOH_ls) <- c("LOHcount","LOHsegments")
  
#HRDetect only needs count of LOH segments
hrd_index <- nrow(LOH_table)

####Structural Variants####
print("summarizing SVs")
SV_bedpe <- try(read.table(SV_bedpe_file,
                           sep = "\t",header = TRUE,
                           stringsAsFactors = FALSE,check.names = FALSE))

if("try-error" %in% class(SV_bedpe) ) {
  
  print("no SVs!")
  SV_ls <- 0
  
} else {
  
  #start sigtools SV analysis
  SV_bedpe_sigtools <- SV_bedpe
  
  #replace svclass name 
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "DEL"] <- "deletion"
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "DUP"] <- "tandem-duplication" 
  SV_bedpe_sigtools$svclass[SV_bedpe_sigtools$svclass == "INV"] <- "inversion" 
  
  #make SV catalog
  SV_catalogue <- bedpeToRearrCatalogue(SV_bedpe_sigtools)
  SV_catalogue_reformat <- sigTools_formatter(input=SV_catalogue,sampleName=sample_name)
  
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
  SV_exp_tissue <- SV_fit_sigs$exposures
  
  #start CHORD SV analysis
  SV_bedpe_CHORD <- SV_bedpe
  SV_bedpe_CHORD$length <-  SV_bedpe_CHORD$start2 - SV_bedpe_CHORD$start1
  
  structuraldf <- SV_bedpe_CHORD[,c("svclass","length")]
  names(structuraldf) <- c("sv_type", "sv_len") 
  #package documentation says col. names aren't necessary, but had trouble running without _these_ names
  
  #extract ALL CHORD signatures
  CHORD_contexts <- extractSigsChord(
    vcf.snv = snvFile_loc,
    vcf.indel = indel_vcf_file,
    df.sv = structuraldf,
    ref.genome=BSgenome.Hsapiens.UCSC.hg38,
    verbose=T
  )
    
  #reformat results for JSON
  CHORD_SV <- as.data.frame(t(CHORD_contexts[1,c(127:145)]))
  unlist(CHORD_contexts[1,c(127:145)])
  
  rearr_catalogue <- as.data.frame(t(SV_catalogue$rearr_catalogue[,1]))
  names(rearr_catalogue) <- rownames(SV_catalogue$rearr_catalogue)
  
  SV_exp_tissue_df <- as.data.frame(t(SV_exp_tissue[1,]))
  
  #make SV list
  SV_ls <- list(nrow(SV_bedpe),CHORD_SV,rearr_catalogue,SV_exp_tissue_df)
  names(SV_ls) <- c("SVcount","SV_CHORD_catalog","SV_sigtools_catalog","SV_sigtools_exposures")

}

####INDEL####
print("summarizing INDELs")

indelTable <- try(read.table(indel_vcf_file,comment.char= "#"))

if("try-error" %in% class(indelTable) | nrow(indelTable) < indelCutoff) {
  
  indel_ls <- 0
  print("no indels!")
  
} else {
  
  #sigtools indel classification
  indels_class <- vcfToIndelsClassification(indel_vcf_file,
                                              sampleID =sample_name,
                                              genome.v = genomeVersion
                                            )
  
  indel_info <- indels_class$count_proportion[c("all.indels","ins","del.mh","del.rep","del.none","del.mh.prop","del.rep.prop","del.none.prop")]
  
  if("try-error" %in% class(SV_bedpe) ) {
    #CHORD will only extract signatures if all files exist
    CHORD_indel <- NA
  } else {
    CHORD_indel <- CHORD_contexts[1,c(97:126)]
  }
  
  #reformat for JSON
  CHORD_indel_df <- as.data.frame(t(CHORD_indel))
  
  #make indel list
  indel_ls <- list(unlist(indel_info[1]),indel_info,CHORD_indel_df)
  names(indel_ls) <- c("indelCount","indel_sigtools_catalog","indel_CHORD_catalog")
  
  #rename indel file for HRDetect input
  names(indel_vcf_file)[1] <- sample_name
  
}

####SNV####
print("summarizing SNVs")

SNVTable <- try(read.table(snvFile_loc,comment.char= "#"))

if("try-error" %in% class(SNVTable) ) {
  
  SNV_ls <- 0
  print("no SNVs!")
  
} else {
  
  #catalog SNVs
  snv_catalogue <- vcfToSNVcatalogue(vcfFilename = snvFile_loc,
                                     genome.v = genomeVersion)
  
  snv_catalogue_reformat <- sigTools_formatter(input=snv_catalogue,
                                               sampleName=sample_name)
  
  #fit COSMIC v2 signatures, aka 'classic' signatures
  print("fitting classic signatures")

  subs_COSMICV2_res <- Fit(catalogues = snv_catalogue_reformat,
                            signatures = COSMIC30_subs_signatures,
                            useBootstrap = TRUE,
                            nboot = boots)
 
  
  #get organ-specific signatures and fit, a la Degasperi et al 2020
  print("fitting tissue signatures")
  
  SNV_sigs <- getOrganSignatures(organ = tissue, 
                                version = "1", 
                                typemut = "subs")
  
  
  SNV_fit_sigs <- Fit(catalogues = snv_catalogue_reformat,
                     signatures = SNV_sigs,
                     useBootstrap = TRUE,
                     nboot = boots)
  
  #add extra step for fitting rare signatures, a la Degasperi et al 2022 
  
  print("fitting rare signatures")
  
  rare_fit <- FitMS(catalogues = snv_catalogue_reformat,"Pancreas",
                         useBootstrap = TRUE,
                         nboot = boots)

  #reformat for JSON
  classic_fit <- as.data.frame(t(subs_COSMICV2_res$exposures[1,]))
  
  tissue_fit <- as.data.frame(t(SNV_fit_sigs$exposures[1,]))
  
  rare_fit_df <-  as.data.frame(t(rare_fit$exposures[1,]))
  
  #make list
  SNV_ls <- list(classic_fit,tissue_fit,rare_fit_df)
  names(SNV_ls) <- c("classic_sigs","tissue_sigs","rare_sigs")
  
}
  
####HRD tests####
print("Performing HRD tests")

if("try-error" %in% class(SNVTable) | "try-error" %in% class(indelTable) |  "try-error" %in% class(SV_bedpe) ) {

  HRD_ls <- NA
  print("data missing, no HRD call!")
  
} else {
  
  #make HRDetect input matrix
  col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  
  input_matrix <- matrix(NA,nrow = 1,
                         ncol = length(col_hrdetect),
                         dimnames = list(sample_name,col_hrdetect))
  
  input_matrix[sample_name,"hrd"] <- hrd_index
  
  #call HRDetect
  HRDetect_res <- HRDetect_pipeline(data_matrix=input_matrix,
                                      bootstrapHRDetectScores=TRUE,
                                      SV_catalogues=SV_catalogue_reformat,
                                      SNV_catalogues=snv_catalogue_reformat,
                                      organ=tissue,
                                      Indels_vcf_files=indel_vcf_file,
                                      genome.v = genomeVersion,
                                      SNV_signature_version = "RefSigv2",
                                      nbootFit=boots
                                      )
  
  #call CHORD
  CHORD_call <- chordPredict(CHORD_contexts, do.bootstrap=T, verbose=T,bootstrap.iters = boots)
  
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
}   

#stick all the results together
all.results <- list(sample_name,LOH_ls,SV_ls,indel_ls,SNV_ls,HR_calls)
names(all.results) <- c("Sample","LOH","SV","indel","SNV","HRD")

#conver to JSON and write
ListJSON <- jsonlite::toJSON(all.results,pretty=TRUE,auto_unbox=TRUE)

write(ListJSON,file = paste(sample_name,".signatures.json",sep=""))


