
####packages####

library(signature.tools.lib)

sigTools_formatter <- function(input,sampleName){
  cat_list <- list()
  names(input)[1] <- "catalogue"
  colnames(input$catalogue) <- sampleName
  cat_list[[1]] <- input$catalogue
  catalogues <- do.call(cbind,cat_list)
  return(catalogues)
}

####start####
#setwd('/Volumes/')

args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
tissue <- args[2]
wkdir <- args[3]

####SNV####
print("reading SNV file... this is the longest step")
snvFile_loc <- paste(wkdir,sample_name,
                     ".filter.deduped.realigned.recalibrated.mutect2.SNP.vcf.gz",sep="")

snv_catalogue <- vcfToSNVcatalogue(vcfFilename = snvFile_loc,genome.v = "hg38")
print("SNV cataloguing done")

snv_catalogue_reformat <- sigTools_formatter(input=snv_catalogue,sampleName=sample_name)

####indel####
indel_vcf_file <- paste(wkdir,sample_name,
                        ".filter.deduped.realigned.recalibrated.mutect2.INDEL.vcf.gz",sep="")
names(indel_vcf_file)[1] <- sample_name

####structural variants####
SV_bedpe <- read.table(paste(wkdir,sample_name,
                               ".somatic.delly.merged.reformat.bedpe",sep=""),
                       sep = "\t",header = TRUE,
                       stringsAsFactors = FALSE,check.names = FALSE)

#replace svclass name (even though documentation suggests this isn't necessary)
SV_bedpe$svclass[SV_bedpe$svclass == "DEL"] <- "deletion"
SV_bedpe$svclass[SV_bedpe$svclass == "DUP"] <- "tandem-duplication" 
SV_bedpe$svclass[SV_bedpe$svclass == "INV"] <- "inversion" 

SV_catalogue <- bedpeToRearrCatalogue(SV_bedpe)

SV_catalogue_reformat <- sigTools_formatter(input=SV_catalogue,sampleName=sample_name)

####LOH#####

ascat.data <- read.table(paste(wkdir,sample_name,"_segments.cna.txt",sep=""),sep="\t",header=TRUE)

hrd_index <- ascatToHRDLOH(ascat.data=ascat.data,SAMPLE.ID=sample_name)

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
                                  #nparallel=threads,
                                  Indels_vcf_files=indel_vcf_file,
                                  genome.v = "hg38"
                                  )

####output####

quantiles <- HRDetect_res$q_5_50_95
names(quantiles) <- c("HRD_low_quant","HRD_median","HRD_top_quant")

print(
  c(
    "HRD_point"=HRDetect_res$hrdetect_output[8],
    quantiles
))

write.table(
  c(
    "HRD_point"=HRDetect_res$hrdetect_output[8],
    quantiles
  ),
  file = paste(wkdir,sample_name,".hrd.txt")
)

#sigsToUse <- getOrganSignatures(tissue,typemut = "subs") 

#subs_fit_res <- SignatureFit_withBootstrap_Analysis(
#  outdir = paste("~/Documents/data/sigtools/signatureFit_snv/",sample_name,
#                 "/",format(Sys.time(), "%Y%m%d%H%M%S"),"/",sep=""),
#  cat = snv_catalogue_reformat,
#  signature_data_matrix = sigsToUse,
#  type_of_mutations = "subs",
#nparallel = threads,
#  nboot = 100

#)

#snv_exp <- convertExposuresFromOrganToRefSigs(expMatrix = subs_fit_res$E_median_filtered,typemut = "subs")

#sigsToUse_rearr <- getOrganSignatures(organ=tissue,typemut = "rearr")

#rearr_fit_res <- 
#  SignatureFit_withBootstrap_Analysis(
#      outdir = paste("~/Documents/data/sigtools/signatureFit_rearr/",sample_name,"/",format(Sys.time(), "%Y%m%d%H%M%S"),"/",sep=""),
#      cat = SV_catalogue_reformat,
#      signature_data_matrix = sigsToUse_rearr,
#      type_of_mutations = "rearr",
#nparallel = threads,
#      nboot = 100

#  )

#SV_exp <- convertExposuresFromOrganToRefSigs(expMatrix = rearr_fit_res$E_median_filtered,typemut = "rearr")

#indel_catalogue <- vcfToIndelsClassification(indelsVCF.file=indel_vcf_file,
#                                             sampleID=sample_name,genome.v = "hg38")

#names(hrd_index) <- "LOH"

#print(
#c(
#  snv_exp[snv_exp != 0,]/sum(snv_exp[snv_exp != 0,]),
#  SV_exp[SV_exp != 0,]/sum(SV_exp[SV_exp != 0,]),
#  hrd_index,
#  "del.mh.prop"=HRDetect_res$indels_classification_table[1,7]
#))


