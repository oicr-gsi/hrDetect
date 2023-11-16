#! /usr/bin/env Rscript

library(signature.tools.lib)
library(BSgenome.Hsapiens.UCSC.hg38)
library(jsonlite)
library(deconstructSigs)
library(dplyr)
library(StructuralVariantAnnotation)

'%ni%' <- function(x,y)!('%in%'(x,y))

call_hrdetect <- function(boots         ,
                          sample_name   ,
                          SV_sigs       ,
                          SNV_sigs      ,
                          snv_df        ,
                          indels_class  ,
                          SV_catalogue  ,
                          seg.data        
){
  
  
  hrdetect_boots <- c("del.mh.prop"  =numeric(),"SNV3"  =numeric(),"SV3"  =numeric(),"SV5"  =numeric(),"hrd"  =numeric(),"SNV8"  =numeric(),"intercept"  =numeric(),
                      "del.mh.prop.w"=numeric(),"SNV3.w"=numeric(),"SV3.w"=numeric(),"SV5.w"=numeric(),"hrd.w"=numeric(),"SNV8.w"=numeric(),"Probability"=numeric())
  
  cat('starting bootstrapping\n0')
  
  for(boot in c(1:boots)){
    if(boot %% 10 == 0){ cat(boot)}else{ cat('.')  }
    
    ## Summarize data
    LOH_ls    <- summarize_LOH(seg.data)
    indel_ls  <- summarize_indels(indels_class)
    SV_ls     <- summarize_SVs(SV_catalogue, SV_sigs, sample_name, nboot = 200)
    SNV_ls    <- summarize_SNVs_dcSigs(snv_df, SNV_sigs, sample_name, subsample_proportion=0.8)
    
    ## Perform test
    input_matrix <- matrix(NA, nrow = 1, ncol = 6,
                            dimnames = list(
                             sample_name,
                             c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
                            )
                           )
    
    input_matrix[sample_name,"del.mh.prop"] <- indel_ls$sim_proportions["del.mh.prop"]
    input_matrix[sample_name,"SNV3"] <- SNV_ls$SBS3
    input_matrix[sample_name,"SNV8"] <- SNV_ls$SBS8
    input_matrix[sample_name,"SV3"]  <- SV_ls$SV_sigtools_exposures$RefSigR3
    input_matrix[sample_name,"SV5"]  <- SV_ls$SV_sigtools_exposures$RefSigR5
    input_matrix[sample_name,"hrd"]  <- unlist(LOH_ls[1])
    
    suppressMessages(
      HRDetect_res <- HRDetect_pipeline(data_matrix = input_matrix, bootstrapHRDetectScores = FALSE)
    )

    boot_results <- c(HRDetect_res$data_matrix, HRDetect_res$hrdetect_output)
    hrdetect_boots <- rbind.data.frame(hrdetect_boots,boot_results )
    names(hrdetect_boots) <- c(colnames(HRDetect_res$data_matrix), paste0(colnames(HRDetect_res$hrdetect_output),".w"))
    
  }   
  cat('done bootstrapping\n')
  
  quantiles <- lapply(hrdetect_boots, function (x) quantile(x, probs = c(5, 50, 95)/100))
  
  return(quantiles)  
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
  
  indel_info <- indels_class$count_proportion[c("all.indels","ins","all.del","del.mh","del.rep","del.none","del.mh.prop","del.rep.prop","del.none.prop")]
  
  sim_proportions <- table(sample(1:3 , indel_info$all.del, replace=TRUE, prob=indel_info[c(7:9)])) / indel_info$all.del
  
  # ensures all three columns exist
  if("1" %ni% names(sim_proportions)){sim_proportions["1"] <- 0}
  if("2" %ni% names(sim_proportions)){sim_proportions["2"] <- 0}
  if("3" %ni% names(sim_proportions)){sim_proportions["3"] <- 0}
  sim_proportions <- sim_proportions[order(names(sim_proportions))]
  
  names(sim_proportions) <-  c("del.mh.prop","del.rep.prop","del.none.prop")
  
  #make indel list
  indel_ls <- list(unlist(indel_info[1]), indel_info, sim_proportions)
  names(indel_ls) <- c("indelCount","indel_sigtools_catalog", "sim_proportions")
  
  return(indel_ls)
}  

summarize_LOH <- function(seg.data, segs_to_sample=40){
  # humans have 40 chromosome arms, so default segs_to_sample is 40
  library(dplyr)
  
  seg.data$length <- seg.data$end - seg.data$start
  sum_of_segment_lengths <- sum(seg.data$length)
  
  LOH.data <- seg.data[seg.data$minorAlleleCopyNumber < 0.5 & seg.data$majorAlleleCopyNumber >= 0.5,]
  sum_of_LOH_lengths <- sum(LOH.data$length)
  
  LOH_probabilities <- c(sum_of_LOH_lengths/sum_of_segment_lengths, 1-sum_of_LOH_lengths/sum_of_segment_lengths)

  LOH_ls <- table(sample(1:2 , segs_to_sample, replace=TRUE, prob=LOH_probabilities))[1] 
  
  names(LOH_ls) <- c("LOHcount")
  
  return(LOH_ls)
  
}

reformat_signatures <- function(SNV_sigs){
  rownames(SNV_sigs) <- SNV_sigs[,1]
  SNV_sigs <- SNV_sigs[,-1]
  SNV_sigs_t <- as.data.frame(t(SNV_sigs))
  return(SNV_sigs_t)
}

summarize_SNVs_dcSigs <- function(snv_df, SNV_sigs, sample_name, subsample_proportion ){
  
  n_snv <- nrow(snv_df)
  
  snv_df <- dplyr::sample_n(snv_df, round(n_snv * subsample_proportion), replace = F)
  
  SNV_sigs_reformatted <- reformat_signatures(SNV_sigs)
  # read input data
  sigs_input <- mut.to.sigs.input(mut.ref=snv_df,
                                  sample.id="Tumor_Sample_Barcode",
                                  chr="Chromosome",
                                  pos="Start_Position",
                                  ref="Reference_Allele",
                                  alt="Allele",
                                  bsg=BSgenome.Hsapiens.UCSC.hg38)
  
  # calculate a sample
  sample_sigs <- whichSignatures(tumor.ref=sigs_input,
                                 signatures.ref=SNV_sigs_reformatted, 
                                 sample.id=sample_name,
                                 contexts.needed=TRUE,
                                 tri.counts.method='exome')
  
  
  # get signature weights
  weights <- sample_sigs$weights
  df_nums <- weights * n_snv
  
  return(df_nums)
  
}

summarize_SNVs_sigTools <- function(snv_catalogue, SNV_sigs, sample_name, giniThresholdScaling=2, threshold_percent=1){
  # Deprecated in .wdl, uses signature.tools.lib version of signatures rather than deconstructSig

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



summarize_SVs <- function(SV_catalogue, SV_sigs, sampleName , nboot = 200, giniThresholdScaling=2, threshold_percent=1){
    

  SV_catalogue_reformat <-  sigTools_formatter(SV_catalogue, sampleName = sample_name)
  
  #fit SV signature
  SV_fit_sigs_all <- Fit(catalogues = SV_catalogue_reformat,
                         signatures = SV_sigs,
                         useBootstrap = TRUE,
                         nboot = nboot,
                         exposureFilterType = "giniScaledThreshold",
                         giniThresholdScaling = giniThresholdScaling, 
                         threshold_percent = threshold_percent
                        )
  
  #extract signature exposures
  SV_exposure <- as.data.frame(t(SV_fit_sigs_all$exposures[1,]))
  
  #make SV list
  SV_ls <- list(SV_catalogue$rearr_catalogue,SV_exposure)
  
  names(SV_ls) <- c("SV_sigtools_catalog","SV_sigtools_exposures")
    
  return(SV_ls)
  
}

simpleEventType <- function(gr) {
  #https://github.com/PapenfussLab/gridss/blob/master/example/simple-event-annotation.R
  pgr = partner(gr)
  
  gr$simple_svtype <-
    ifelse(seqnames(gr) != seqnames(pgr), "translocation", # inter-chromosomosal
         ifelse(strand(gr) == strand(pgr), "inversion",
                ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
                       ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "deletion",
                              "tandem-duplication"))))
  return(gr)
}

SV_vcf_cataloger <- function(SV_vcf_location, sample_name, genome="hg38"){
  
  vcf <- VariantAnnotation::readVcf(SV_vcf_location, genome)

  gr <- StructuralVariantAnnotation::breakpointRanges(vcf, nominalPosition=T, inferMissingBreakends=TRUE)
  
  gr <- simpleEventType(gr)
  
  simplebed <- data.frame(
    chrom1 = seqnames(gr),
    start1 = start(gr),
    end1 = end(gr),
    chrom2 = seqnames(partner(gr)),
    start2 = start(partner(gr)),
    end2 = end(partner(gr)),
    sample = sample_name,
    svclass = gr$simple_svtype,
    score = gr$FILTER,
    strand = "."
  )
  
  simplebed <- simplebed %>% filter(svclass != "translocation")
  # Just the lower of the two breakends, so we don't output everything twice
  simplebed <- simplebed[order(simplebed$chrom1, simplebed$start1, simplebed$chrom2, simplebed$start2),]
  simplebed <- simplebed[!duplicated(simplebed),]

  SV_catalogue <- bedpeToRearrCatalogue2(sv_bedpe=simplebed)
  return(SV_catalogue)
} 

bedpeToRearrCatalogue2 <- function(sv_bedpe, kmin = 10, PEAK.FACTOR = 10) {
  #Replaces bedpeToRearrCatalogue from signatures.tools which crashed if it filtered out all the SVs
  library(signature.tools.lib)
  # this is used for per-sample clustering of both single-base substitutions and rearrangement breakpoints
  
  rearrangement.clustering_bedpe <- function(sv_bedpe,
                                             plot.path=NA,
                                             kmin=10,# how many points at minimum in a peak, for the pcf algorithm
                                             kmin.samples=kmin, # how many different samples at minimum in  a peak
                                             gamma.sdev=25, #
                                             PEAK.FACTOR=4,
                                             thresh.dist=NA,
                                             gamma=NA,
                                             kmin.filter=kmin # if the pcf parameter is different from the definition of a peak
  ) {
    
    #add an id to the rearrangement
    sv_bedpe$id <- 1:nrow(sv_bedpe)
    
    #functions below expect rows to be organised by chromosomes and ordered by position on the chromosome
    
    #prepare a dataframe for the calculation
    rearrs.left <- sv_bedpe[,c('chrom1','start1','sample')]
    names(rearrs.left ) <- NA
    rearrs.right <- sv_bedpe[,c('chrom2','start2','sample')]
    names(rearrs.right ) <- NA
    rearrs.cncd <- rbind(rearrs.left , rearrs.right  )
    colnames(rearrs.cncd) <- c('chr', 'position', 'sample')
    rearrs.cncd$isLeft <- c(rep(TRUE, nrow(rearrs.left)), rep(FALSE, nrow(rearrs.left)))
    rearrs.cncd$id <- c(sv_bedpe$id, sv_bedpe$id)
    # sample.bps <- rearrs.cncd
    #need to reorder
    sample.bps <- NULL
    for (chrom_i in unique(rearrs.cncd$chr)){
      tmptab <- rearrs.cncd[rearrs.cncd$chr==chrom_i,,drop=FALSE]
      tmptab <- tmptab[order(tmptab$position),,drop=FALSE]
      sample.bps <- rbind(sample.bps,tmptab)
    }
    rownames(sample.bps) <- 1:nrow(sample.bps)
    
    #run the algorithm
    genome.size <- 3 * 10^9
    MIN.BPS <- 10 # minimal number of breakpoints on a chromosome to do any any segmentation
    
    logScale <- FALSE
    
    exp.dist <-genome.size/nrow(sample.bps)
    
    if (logScale) {
      sample.bps$intermut.dist <- log10(calcIntermutDist(sample.bps, first.chrom.na=FALSE)$distPrev) # calculate the distances between the breakpoints
      if (is.na(thresh.dist)) {
        thresh.dist <- log10(exp.dist/PEAK.FACTOR) # calculate the threshold to call a peak
      }
    } else {
      
      sample.bps$intermut.dist <- calcIntermutDist(sample.bps, first.chrom.na=FALSE)$distPrev
      if (is.na(thresh.dist)) {
        thresh.dist <- exp.dist/PEAK.FACTOR
      }
    }
    
    
    if (is.na(gamma) & !is.na(gamma.sdev)) {
      # compute the mean absolute deviation
      sdev <- getMad(sample.bps$intermut.dist);
      gamma <- gamma.sdev*sdev
    }
    
    
    
    sample.bps$is.clustered.single <- rep(FALSE, nrow(sample.bps))
    
    all.kat.regions <- data.frame()
    
    for (chrom in unique(sample.bps$chr)) { # loop over chromosomes
      
      sample.bps.flag <- sample.bps$chr==chrom #   breakpoints on a current chromosome
      #
      if (sum(sample.bps.flag )>MIN.BPS ) { # if there are enough breakpoints on a chromosome to run pcf
        
        data.points <- sample.bps$intermut.dist[sample.bps.flag]
        
        res = exactPcf(data.points, kmin, gamma, T)
        
        #reorder results
        sample.bps$mean.intermut.dist[sample.bps.flag] <- res$yhat
        
        # prepare the points for pcf
        subs <- data.frame(chr=sample.bps$chr[sample.bps.flag], pos=sample.bps$position[sample.bps.flag], sample=sample.bps$sample[sample.bps.flag])
        kat.regions <- extract.kat.regions(res, thresh.dist, subs, doMerging=TRUE, kmin.samples=1,  kmin.filter= kmin.filter) # extract peaks, this is special case as we want at least kmin samples
        
        all.kat.regions <- rbind(all.kat.regions, kat.regions)
        if (!is.null(kat.regions) && nrow( kat.regions )>0) { # if there are any kataegis regions found on this chormosome
          for (k in 1:nrow(kat.regions)) {
            
            sample.bps$is.clustered.single[which(sample.bps.flag)[ kat.regions$firstBp[k] : kat.regions$lastBp[k]]] <- TRUE # call all breakpoints as clustered
          }
        }
      } else {
        
        sample.bps$mean.intermut.dist[sample.bps.flag] <- mean(sample.bps$intermut.dist[sample.bps.flag])
      }
    }
    
    
    
    if (!logScale) { # even if pcf was run on non-logged distances, I log the output
      sample.bps$intermut.dist <- log10(sample.bps$intermut.dist)
      sample.bps$mean.intermut.dist <- log10(sample.bps$mean.intermut.dist)
    }
    
    # a rearrangement is in a cluster if any of its breakpoints are
    sample.bps$is.clustered <- sample.bps$is.clustered.single
    sample.bps$is.clustered[sample.bps$id %in% subset(sample.bps, is.clustered.single==TRUE)$id] <- TRUE
    
    # mark both breakpoints of a rearrangement as clustered if any is
    sv_bedpe$is.clustered <- sv_bedpe$id %in% sample.bps$id[sample.bps$is.clustered]
    
    result <- list()
    result$sv_bedpe <- sv_bedpe
    result$clustering_regions <- all.kat.regions
    result
  }
  
  
  prepare.rearr.catalogue_fromAnnotatedBedpe <- function(sv_bedpe) {
    
    catalogue.labels <- c('clustered_del_1-10Kb', 'clustered_del_10-100Kb', 'clustered_del_100Kb-1Mb', 'clustered_del_1Mb-10Mb', 'clustered_del_>10Mb', 'clustered_tds_1-10Kb', 'clustered_tds_10-100Kb', 'clustered_tds_100Kb-1Mb', 'clustered_tds_1Mb-10Mb', 'clustered_tds_>10Mb', 'clustered_inv_1-10Kb', 'clustered_inv_10-100Kb', 'clustered_inv_100Kb-1Mb', 'clustered_inv_1Mb-10Mb', 'clustered_inv_>10Mb', 'clustered_trans', 'non-clustered_del_1-10Kb', 'non-clustered_del_10-100Kb', 'non-clustered_del_100Kb-1Mb', 'non-clustered_del_1Mb-10Mb', 'non-clustered_del_>10Mb', 'non-clustered_tds_1-10Kb', 'non-clustered_tds_10-100Kb', 'non-clustered_tds_100Kb-1Mb', 'non-clustered_tds_1Mb-10Mb', 'non-clustered_tds_>10Mb', 'non-clustered_inv_1-10Kb', 'non-clustered_inv_10-100Kb', 'non-clustered_inv_100Kb-1Mb', 'non-clustered_inv_1Mb-10Mb', 'non-clustered_inv_>10Mb', 'non-clustered_trans')
    
    all_catalogues <- as.data.frame(matrix(nrow = length(catalogue.labels),ncol = 0))
    rownames(all_catalogues) <- catalogue.labels
    
    updated_sv_bedpe <- NULL
    
    if (nrow(sv_bedpe)>0){
      for (sample_name in unique(sv_bedpe$sample)){
        sample.rearrs <- sv_bedpe[sv_bedpe$sample==sample_name,]
        
        rearr_catalogue <- as.data.frame(matrix(0,nrow = length(catalogue.labels),ncol = 1))
        
        if (nrow(sample.rearrs)>0) {
          
          label1 <- rep('non-clustered', nrow(sample.rearrs))
          label1[ sample.rearrs$is.clustered] <- 'clustered'
          
          label2 <- rep('', nrow(sample.rearrs))
          label2[ sample.rearrs$svclass=='deletion'] <- '_del'
          label2[ sample.rearrs$svclass=='translocation'] <- '_trans'
          label2[ sample.rearrs$svclass=='inversion'] <- '_inv'
          label2[ sample.rearrs$svclass=='tandem-duplication'] <- '_tds'
          
          label3 <- rep('', nrow(sample.rearrs))
          sample.rearrs$bkdist <- abs(sample.rearrs$start2 - sample.rearrs$start1)
          label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist<=1e4] <- '_1-10Kb'
          label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e4 & sample.rearrs$bkdist<=1e5 ] <- '_10-100Kb'
          label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e5 & sample.rearrs$bkdist<=1e6 ] <- '_100Kb-1Mb'
          label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e6 & sample.rearrs$bkdist<=1e7 ] <- '_1Mb-10Mb'
          label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e7 ] <- '_>10Mb'
          
          sample.rearrs$catalogue.label <- paste0(label1, label2, label3)
          updated_sv_bedpe <- rbind(updated_sv_bedpe,sample.rearrs)
          
          sample.table <- as.data.frame(table(sample.rearrs$catalogue.label),drop=FALSE)
          rownames(sample.table) <- sample.table$Var
          
          rearr_catalogue <-  sample.table[as.character(catalogue.labels), 'Freq',drop=FALSE ]
          
        }
        
        rearr.catalogue <- rearr_catalogue
        rownames(rearr.catalogue) <- catalogue.labels
        colnames(rearr.catalogue) <- sample_name
        rearr.catalogue[is.na(rearr.catalogue)] <- 0
        
        all_catalogues <- cbind(all_catalogues,rearr.catalogue)
      }
    }else{
      all_catalogues <- as.data.frame(matrix(0,nrow = length(catalogue.labels),ncol = 1))
      rownames(all_catalogues) <- catalogue.labels
      updated_sv_bedpe <- sv_bedpe
    }
    
    returnObj <- list()
    returnObj$SV_catalogue <- all_catalogues
    returnObj$updated_sv_bedpe <- updated_sv_bedpe
    return(returnObj)
    
  }
  
  
  
  classifyRearrangementsFromBedpe <- function(sv_bedpe){
    svclass <- c()
    for (i in 1:nrow(sv_bedpe)){
      if(sv_bedpe[i,"chrom1"]!=sv_bedpe[i,"chrom2"]){
        svclass <- c(svclass,"translocation")
      }else if(sv_bedpe[i,"strand1"]!=sv_bedpe[i,"strand2"]){
        svclass <- c(svclass,"inversion")
      }else if(sv_bedpe[i,"strand1"]=="+"){
        svclass <- c(svclass,"deletion")
      }else if(sv_bedpe[i,"strand1"]=="-"){
        svclass <- c(svclass,"tandem-duplication")
      }
    }
    sv_bedpe[,"svclass"] <- svclass
    #return updated df
    sv_bedpe
  }
  
  
  build_junctions_catalogue <- function(annotated_bedpe){
    
    sizes <- c("_1-3",
               "_4-10",
               "_>10")
    junctions_catalogue_channels <- paste(rep(c("clustered","non-clustered"),each=7),
                                          rep(c(rep("_non-templated",3),"_other",rep("_homologous",3)),2),
                                          rep(c(rev(sizes),"",sizes),2),sep = "")
    junctions_catalogue <- data.frame(sample=rep(0,length(junctions_catalogue_channels)),
                                      row.names = junctions_catalogue_channels,
                                      stringsAsFactors = F)
    if(nrow(annotated_bedpe)==0){
      # cut short here returning a 0 catalogue if there are no SVs
      return(junctions_catalogue)
    }
    
    # # check if non-template and micro-homology columns are present
    if(all(c("non-template","micro-homology") %in% colnames(annotated_bedpe))){
      for (clustered in c(TRUE,FALSE)){
        # clustered <- F
        channelc <- ifelse(clustered,"clustered","non-clustered")
        for (typebp in c("_non-templated","_homologous","_other")){
          # typebp <- "_non-templated"
          channel <- paste0(channelc,typebp)
          if(typebp %in% c("_non-templated","_homologous")){
            bedpecol <- ifelse(typebp=="_non-templated","non-template","micro-homology")
            currentseqlength <- nchar(annotated_bedpe[annotated_bedpe[,bedpecol]!="." & annotated_bedpe$is.clustered==clustered,bedpecol])
            if(length(currentseqlength)>0){
              currentlengthtable <- table(currentseqlength)
              for(i in 1:length(currentlengthtable)){
                # i <- 1
                tmpn <- as.numeric(names(currentlengthtable))[i]
                if(tmpn>0 & tmpn<=3){
                  junctions_catalogue[paste0(channel,"_1-3"),"sample"] <- junctions_catalogue[paste0(channel,"_1-3"),"sample"] + currentlengthtable[i]
                }else if(tmpn>3 & tmpn<=10){
                  junctions_catalogue[paste0(channel,"_4-10"),"sample"] <- junctions_catalogue[paste0(channel,"_4-10"),"sample"] + currentlengthtable[i]
                }else if(tmpn>10){
                  junctions_catalogue[paste0(channel,"_>10"),"sample"] <- junctions_catalogue[paste0(channel,"_>10"),"sample"] + currentlengthtable[i]
                }
              }
            }
          }else{
            tmpn <- sum(annotated_bedpe[,"non-template"]=="." &  annotated_bedpe[,"micro-homology"]=="." & annotated_bedpe$is.clustered==clustered)
            junctions_catalogue[channel,"sample"] <- junctions_catalogue[channel,"sample"] + tmpn
          }
        }
      }
      return(junctions_catalogue)
    }else{
      # return NULL if there are SVs but no "non-template" and "micro-homology" columns
      return(NULL)
    }
    
  }
  
  required_cols <- c("chrom1", "start1", "end1", "chrom2", 
                     "start2", "end2", "sample")
  if (!length(intersect(required_cols, colnames(sv_bedpe))) == 
      length(required_cols) & ("svclass" %in% colnames(sv_bedpe) | 
                               all(c("strand1", "strand2") %in% colnames(sv_bedpe)))) {
    stop("[error bedpeToRearrCatalogue] missing columns in subs data frame, following columns required: chrom1, start1, end1, chrom2, start2, end2 and sample. In addition either svclass or strand1 and strand2. Check ?bedpeToRearrCatalogue for details.")
  }
  clusters_table <- NULL
  if (nrow(sv_bedpe) > 0) {
    if (length(unique(sv_bedpe$sample)) > 1) {
      tmpsamplenames <- unique(sv_bedpe$sample)
      tmpcounts <- table(sv_bedpe$sample)
      pos <- which.max(tmpcounts)
      samplechoice <- names(pos)
      sv_bedpe <- sv_bedpe[sv_bedpe$sample == samplechoice, 
                           , drop = F]
      message("[warning bedpeToRearrCatalogue] bedpeToRearrCatalogue sample column should have only one sample name, however multiple sample names were detected: ", 
              paste(tmpsamplenames, collapse = ", "), ". Trying to fix this by using only the sample with the largest number of rearrangements (", 
              samplechoice, "). ", "This fix should work if there are a few rearrangements from the germline sample, so all we would do is to remove the germline sample name and variants ", 
              "while keeping all the tumour sample variants.")
    }
  }
  if (nrow(sv_bedpe) > 0) {
    if (!"svclass" %in% colnames(sv_bedpe)) {
      if ("strand1" %in% colnames(sv_bedpe) & "strand2" %in% 
          colnames(sv_bedpe)) {
        sv_bedpe <- classifyRearrangementsFromBedpe(sv_bedpe)
      }
      else {
        message("[error bedpeToRearrCatalogue] cannot classify rearrangements: svclass column missing, and cannot compute it because strand1 and strand2 are missing.")
        return(NULL)
      }
    }
    all_sv_annotated <- sv_bedpe
    bkdist <- abs(sv_bedpe$start2 - sv_bedpe$start1)
    sv_bedpe[sv_bedpe$svclass != "translocation", "length"] <- bkdist[sv_bedpe$svclass != 
                                                                        "translocation"]
    toberemoved <- sv_bedpe$svclass != "translocation" & 
      bkdist < 1000
      
    if (sum(toberemoved) > 0) {
      message("[warning bedpeToRearrCatalogue] ignoring rearrangements shorter than 1kb (", 
              sum(toberemoved), " out of ", nrow(all_sv_annotated), 
              ")")
      all_sv_annotated[toberemoved, "FILTER"] <- "length<1e3"
      all_sv_annotated[!toberemoved, "FILTER"] <- "PASS"
      sv_bedpe <- sv_bedpe[!toberemoved, ]
    }
    if (sum(toberemoved) != nrow(all_sv_annotated)) {
      
      clustering.result <- rearrangement.clustering_bedpe(sv_bedpe, 
                                                          plot.path = NA, kmin = kmin, kmin.samples = 1, gamma.sdev = 25, 
                                                          PEAK.FACTOR = PEAK.FACTOR, thresh.dist = NA)
      sv_bedpe <- clustering.result$sv_bedpe
      clustering_regions <- clustering.result$clustering_regions
    } else {
      message("[warning bedpeToRearrCatalogue] ignoring rearrangements shorter than 1kb (", 
              sum(toberemoved), " out of ", nrow(all_sv_annotated), 
              ") -- leaving none")
      clustering_regions <- NULL
    }
  }
  else {
    all_sv_annotated <- NULL
    clustering_regions <- NULL
  }
  resSVcat <- prepare.rearr.catalogue_fromAnnotatedBedpe(sv_bedpe)
  rearr_catalogue <- resSVcat$SV_catalogue
  sv_bedpe <- resSVcat$updated_sv_bedpe
  junctions_catalogue <- build_junctions_catalogue(sv_bedpe)
  if (!is.null(junctions_catalogue)) 
    colnames(junctions_catalogue) <- colnames(rearr_catalogue)
  if (nrow(sv_bedpe) > 0) {
    annotated_bedpe <- sv_bedpe
  }
  else {
    annotated_bedpe <- NULL
  }
  return_list <- list()
  return_list$rearr_catalogue <- rearr_catalogue
  return_list$junctions_catalogue <- junctions_catalogue
  return_list$annotated_bedpe <- annotated_bedpe
  return_list$all_sv_annotated <- all_sv_annotated
  if (!is.null(clustering_regions)) {
    if (nrow(clustering_regions) > 0) 
      return_list$clustering_regions <- clustering_regions
  }
  return(return_list)
}



