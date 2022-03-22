
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

#setwd('/Volumes/')

args = commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
tissue <- args[2]
snvFile_loc  <- args[3]
indel_vcf_file <- args[4]
SV_bedpe_file <- args[5]
LOH_seg_file <- args[6]

####Import files####

##SNV##


snv_catalogue <- vcfToSNVcatalogue(vcfFilename = snvFile_loc,genome.v = "hg38")

snv_catalogue_reformat <- sigTools_formatter(input=snv_catalogue,sampleName=sample_name)

##indel##


names(indel_vcf_file)[1] <- sample_name

##structural variants##

SV_bedpe <- read.table(SV_bedpe_file,
                       sep = "\t",header = TRUE,
                       stringsAsFactors = FALSE,check.names = FALSE)

#replace svclass name (even though documentation suggests this isn't necessary)
SV_bedpe$svclass[SV_bedpe$svclass == "DEL"] <- "deletion"
SV_bedpe$svclass[SV_bedpe$svclass == "DUP"] <- "tandem-duplication" 
SV_bedpe$svclass[SV_bedpe$svclass == "INV"] <- "inversion" 

SV_catalogue <- bedpeToRearrCatalogue(SV_bedpe)

SV_catalogue_reformat <- sigTools_formatter(input=SV_catalogue,sampleName=sample_name)

##LOH##

ascat.data <- read.table(LOH_seg_file,sep="\t",header=TRUE)
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
                                  Indels_vcf_files=indel_vcf_file,
                                  genome.v = "hg38"
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



