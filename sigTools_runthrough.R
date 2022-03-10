
####packages####
##test install##
#setwd('~/Documents/GitHub/signature.tools.lib')
#devtools::test()

library(signature.tools.lib)

####start####
setwd('/Volumes/')

#args = commandArgs(trailingOnly=TRUE)
#sample_name <- args[1]

#sample_name1 <- "OCT_010638_Ov_P"
#sample_name2 <- ""
#tissue="Ovary"

sample_name1 <- "PANX_1312_Lv_M"
sample_name2 <- "_100-PM-035_LCM6"
tissue="Liver"

sample_name <- paste(sample_name1,sample_name2,sep="")
 

####SNV####
snvFile_loc <- paste("cgi/scratch/fbeaudry/",sample_name1,"_WG",sample_name2,".filter.deduped.realigned.recalibrated.mutect2.SNP.vcf.gz",sep="")

snv_catalogue <- vcfToSNVcatalogue(vcfFilename = snvFile_loc,genome.v = "hg38")

SNVcat_list <- list()
colnames(snv_catalogue$catalogue) <- sample_name
SNVcat_list[[1]] <- snv_catalogue$catalogue
SNV_catalogues <- do.call(cbind,SNVcat_list)

sigsToUse <- getOrganSignatures(tissue,typemut = "subs") 

subs_fit_res <- SignatureFit_withBootstrap_Analysis(
  outdir = paste("~/Documents/data/sigtools/signatureFit_snv/",sample_name,"/",format(Sys.time(), "%Y%m%d%H%M%S"),"/",sep=""),
  cat = SNV_catalogues,
  signature_data_matrix = sigsToUse,
  type_of_mutations = "subs",
  nboot = 100,
  nparallel = 2
)

#convert tissue to COSMIC
snv_exp <- convertExposuresFromOrganToRefSigs(expMatrix = subs_fit_res$E_median_filtered,typemut = "subs")

####structural variants####
SV_bedpe_file <- c(paste("cgi/scratch/fbeaudry/",sample_name1,"_WG",sample_name2,".somatic.delly.merged.reformat.bedpe",sep=""))
#names(SV_bedpe_file)[1] <- sample_name

SV_bedpe <- read.table(SV_bedpe_file,sep = "\t",header = TRUE,
                       stringsAsFactors = FALSE,check.names = FALSE)

SV_bedpe$svclass[SV_bedpe$svclass == "DEL"] <- "deletion"
SV_bedpe$svclass[SV_bedpe$svclass == "DUP"] <- "tandem-duplication" #check that every DUP is tandem
SV_bedpe$svclass[SV_bedpe$svclass == "INV"] <- "inversion" #check that every DUP is tandem


SV_catalogue <- bedpeToRearrCatalogue(SV_bedpe)
SVcat_list <- list()
colnames(SV_catalogue$rearr_catalogue) <- sample_name
SVcat_list[[1]] <- SV_catalogue$rearr_catalogue
SV_catalogues <- do.call(cbind,SVcat_list)

sigsToUse_rearr <- getOrganSignatures(organ=tissue,typemut = "rearr")

   
rearr_fit_res <- SignatureFit_withBootstrap_Analysis(
                                                    outdir = paste("~/Documents/data/sigtools/signatureFit_rearr/",sample_name,"/",format(Sys.time(), "%Y%m%d%H%M%S"),"/",sep=""),
                                                    cat = SV_catalogues,
                                                    signature_data_matrix = sigsToUse_rearr,
                                                    type_of_mutations = "rearr",
                                                    nboot = 100,
                                                    nparallel = 2
                                                    )

sv_exp <- convertExposuresFromOrganToRefSigs(expMatrix = rearr_fit_res$E_median_filtered,typemut = "rearr")

#rearrangement_table <- fread("/Library/Frameworks/R.framework/Resources/library/signature.tools.lib/data/2019_01_10_ConversionMatrix_rearr.tsv")


####LOH#####

CNV_tab_file <- paste("cgi/scratch/fbeaudry/",sample_name1,"_WG",sample_name2,"_segments.cna.txt",sep="")
names(CNV_tab_file)[1] <- sample_name

ascat.data <- read.table(CNV_tab_file,sep="\t",header=TRUE)

hrd_index <- ascatToHRDLOH(ascat.data=ascat.data,SAMPLE.ID=sample_name)

####indel####
indel_vcf_file <- paste("cgi/scratch/fbeaudry/",sample_name1,"_WG",sample_name2,".filter.deduped.realigned.recalibrated.mutect2.indels.vcf.gz",sep="")
names(indel_vcf_file)[1] <- sample_name
indel_catalogue <- vcfToIndelsClassification(indelsVCF.file=indel_vcf_file,sampleID=sample_name,genome.v = "hg38")

####HRD test####
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = 1,
                       ncol = length(col_hrdetect),
                       dimnames = list(sample_name,col_hrdetect))

input_matrix[sample_name,"SNV3"] <- snv_exp["Ref.Sig.3",1]
input_matrix[sample_name,"SNV8"] <- snv_exp["Ref.Sig.8",1]

input_matrix[sample_name,"SV3"] <- sv_exp["RefSig R3",1]

input_matrix[sample_name,"SV5"] <- sv_exp["RefSig R5",1]
#input_matrix[sample_name,"SV5"] <- sv_exp["RefSig R9",1] #error in pancreas database??

input_matrix[sample_name,"hrd"] <- hrd_index

input_matrix[sample_name,"del.mh.prop"] <- indel_catalogue$count_proportion[1,7]

input_matrix

HRDetect_res <- HRDetect_pipeline(data_matrix=input_matrix, nparallel = 2)

HRDetect_res$hrdetect_output[8]





