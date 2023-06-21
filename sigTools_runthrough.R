#! /usr/bin/env Rscript

library(optparse)

option_list = list(
  make_option(c("-s", "--sampleName"), type="character", default=NULL, help="sample name", metavar="character"),
  make_option(c("-r", "--SVrefSigs"), type="character", default=NULL, help="oncotree code", metavar="character"),
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
SVrefSigs           <-  opt$SVrefSigs
snv_vcf_location    <-  opt$snvFile
indel_vcf_location  <-  opt$indelFile
SV_bedpe_location   <-  opt$SVFile
LOH_seg_location    <-  opt$LOHFile

call_hrdetect(boots               ,
                          genomeVersion       ,
                          indelCutoff         ,
                          snvCutoff           ,
                          sample_name         ,
                          SVrefSigs           ,
                          snv_vcf_location    ,
                          indel_vcf_location  ,
                          SV_bedpe_location   ,
                          LOH_seg_location    
)
  
