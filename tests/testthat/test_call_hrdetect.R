#! /usr/bin/env Rscript

library(testthat)

basedir <- paste(Sys.getenv(c("BASE_DIR")), sep='/')
source(paste0(basedir, "/call_hrdetect.R"))

sample_name <- "TUMOUR_SAMPLE"
indel_vcf_location <- paste0(basedir, "tests/testthat/data/indel.VAF.vcf.gz")
indel_ONE_vcf_location <- paste0(basedir, "tests/testthat/data/indel.ONE.vcf.gz")

snv_df <- read.table(paste0(basedir, "tests/testthat/data/snp.VAF.vcf.gz"), sep="\t", header=TRUE)
names(snv_df)[c(1,2,4,5)] <- c("Chromosome", "Start_Position", "Reference_Allele", "Allele")
snv_df$Tumor_Sample_Barcode <- sample_name

SV_sigs  <- read.table(paste0(basedir,"data/RefSigv0_Rearr.tsv"), sep = "\t", header = TRUE)
SNV_sigs <- read.table(paste0(basedir,"data/COSMIC_v1_SBS_GRCh38.txt"), sep = "\t", header = TRUE)

SV_vcf_location <- paste0(basedir, "tests/testthat/data/structural.PASS.vcf")
SV_ONE_vcf_location <- paste0(basedir, "tests/testthat/data/structural.ONE.vcf")


seg.data <- read.table(paste0(basedir,"tests/testthat/data/cnv.somatic.tsv"), sep="\t", header=TRUE)

##

test_that("call_hrdetect can handle a bunch of nulls",
  {
    
    quantiles <-
      call_hrdetect(boots=10       ,
                    sample_name   ,
                    SV_sigs       ,
                    SNV_sigs      ,
                    NULL        ,
                    NULL  ,
                    NULL  ,
                    seg.data        
      )
    expect_equal(round(quantiles$Probability.w[[2]],2),0)
  }
)

test_that("hrdetect can handle a bunch of zeroes",
          {
            input_matrix <- matrix(c(0,0,0,0,0,0), nrow = 1, ncol = 6,
                                   dimnames = list(
                                       sample_name,
                                       c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
                                     )
                                  
                                   )
            HRDetect_res <- HRDetect_pipeline(data_matrix = input_matrix, bootstrapHRDetectScores = FALSE)
            expect_equal(round(HRDetect_res$hrdetect_output[[8]],2),0)
          }
)


test_that("summarize_indels can handle cases with 1 indel",
          {
            suppressWarnings(
              indels_class    <- vcfToIndelsClassification(indel_ONE_vcf_location, sample_name, genome.v = "hg38")
            )
            indel_ls <- summarize_indels(indels_class, seed=12345)
            del.mh.prop.obs <- indel_ls$sim_proportions["del.mh.prop"]
            expect_equal(del.mh.prop.obs[[1]], 0)
            
          }
)

test_that("summarize_SVs can summarize structural variants",
          {
            suppressWarnings(
              SV_catalogue  <- SV_vcf_cataloger(SV_vcf_location, sample_name)
            )
            SV_ls     <- summarize_SVs(SV_catalogue, SV_sigs, sample_name, nboot = 200, seed=12345)

            expect_equal(SV_ls$SV_sigtools_exposures$RefSigR3, 0)
            expect_equal(round(SV_ls$SV_sigtools_exposures$RefSigR5, 2), 37.08)
            
          }
)

test_that("summarize_SVs can handle cases with 1 SV",
          {
            
            SV_catalogue  <- SV_vcf_cataloger(SV_ONE_vcf_location, sample_name)
            
            SV_ls     <- summarize_SVs(SV_catalogue, SV_sigs, sample_name, nboot = 200, seed=12345)
            
            expect_equal(SV_ls$SV_sigtools_exposures$RefSigR5, 2)
            
          }
)


test_that("summarize_indels summarizes indels",
          {
            suppressWarnings(
              indels_class    <- vcfToIndelsClassification(indel_vcf_location, sample_name, genome.v = "hg38")
            )
            indel_ls <- summarize_indels(indels_class, seed=12345)
            del.mh.prop.obs <- indel_ls$sim_proportions["del.mh.prop"]
            del.mh.prop.exp <- c("del.mh.prop" = 0.3553719)
            expect_equal(del.mh.prop.obs, del.mh.prop.exp)
            
          }
)



test_that("summarize_SNVs_dcSigs can summarize and deconstruct signatures",
          {
            SNV_ls    <- summarize_SNVs_dcSigs(snv_df, SNV_sigs, sample_name, subsample_proportion=0.8, seed=12345)
            expect_equal( round(SNV_ls$SBS8,2), 70.13)
            expect_equal( round(SNV_ls$SBS3,2), 0)
            
          }
)

test_that("summarize_SNVs_dcSigs can handle cases with 1 SNVs",
          {
            
            captured_warning <- capture_warning(SNV_ls    <- summarize_SNVs_dcSigs(snv_df[1,], SNV_sigs, sample_name, subsample_proportion=0.8, seed=12345))
            expected_warning <- list("message"= "Some samples have fewer than 50 mutations:\n  TUMOUR_SAMPLE")
            expect_equal(round(SNV_ls$SBS3,2), 0)
            expect_equal(captured_warning[1], expected_warning)
          }
)
