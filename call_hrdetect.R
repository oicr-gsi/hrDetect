#! /usr/bin/env Rscript

library(signature.tools.lib)
library(BSgenome.Hsapiens.UCSC.hg38)
library(jsonlite)
library(deconstructSigs)
library(dplyr)
library(StructuralVariantAnnotation)

'%ni%' <- function(x,y)!('%in%'(x,y))

bedpeToRearrCatalogue2 <- function(sv_bedpe, kmin = 10, PEAK.FACTOR = 10) {
  #Replaces bedpeToRearrCatalogue from signatures.tools which crashed if it filtered out all the SVs
  require(signature.tools.lib)
  # this is used for per-sample clustering of both single-base substitutions and rearrangement breakpoints
  
  assignPvalues <- function(kat.regions, chrom.bps, bp.rate=NA) {
    
    if (is.na(bp.rate)) { # estimate the chromosome rate
      left.bp <- min(chrom.bps$pos)
      right.bp <-  max(chrom.bps$pos)
      bp.rate <- nrow(chrom.bps)/ (right.bp - left.bp)
    }
    
    # assume binomial distribution
    kat.regions$pvalue <- 1-pbinom(kat.regions$number.bps, kat.regions$end.bp - kat.regions$start.bp, bp.rate)
    
    kat.regions$d.seg<- (kat.regions$number.bps/( kat.regions$end.bp - kat.regions$start.bp))
    
    kat.regions$rate.factor <- kat.regions$d.seg/bp.rate
    
    kat.regions
  }
  
  hotspotInfo <- function(kat.regions.all, subs, segInterDist=c()) {
    if(nrow(kat.regions.all)>0){
      for(r in 1:nrow(kat.regions.all)){
        
        # indices of the breakpoints in the hotspot
        subs.hotspot <-subs[kat.regions.all$firstBp[r]:kat.regions.all$lastBp[r],]
        
        kat.regions.all[r,'start.bp'] <- min(subs.hotspot$pos)
        kat.regions.all[r,'end.bp'] <- max(subs.hotspot$pos)
        kat.regions.all[r,'length.bp'] <-  kat.regions.all[r,'end.bp'] - kat.regions.all[r,'start.bp'] 
        kat.regions.all[r,'number.bps'] <- nrow(subs.hotspot)
        kat.regions.all[r,'number.bps.clustered'] <- sum(subs.hotspot$is.clustered)
        
        if (length(segInterDist)>0 & is.na(kat.regions.all[r,'avgDist.bp'])) {
          kat.regions.all[r,'avgDist.bp'] <- mean(segInterDist[kat.regions.all$firstBp[r]:kat.regions.all$lastBp[r]])
        }
        kat.regions.all[r,'no.samples'] <- length(unique(subs.hotspot$sample))
        
        if ('pf' %in% colnames(subs.hotspot)){
          kat.regions.all[r,'no.del'] <- nrow(subset(subs.hotspot, pf==2))
          kat.regions.all[r,'no.dup'] <- nrow(subset(subs.hotspot, pf==4))
          kat.regions.all[r,'no.inv'] <- nrow(subset(subs.hotspot, pf==1 | pf==8))
          kat.regions.all[r,'no.trn'] <- nrow(subset(subs.hotspot, pf==32))
        }
        
      } # for all peaks
    } # if there is at least one peak
    kat.regions.all
  }
  
  
  extract.kat.regions <- function (res, imd, subs,  kmin.samples=10, pvalue.thresh=1, rate.factor.thresh=1, doMerging=FALSE, kmin.filter=NA, bp.rate=NA) {
    
    segInterDist <-  res$yhat
    kataegis.threshold <- imd
    
    kat.regions.all = data.frame()	
    
    chr <- as.character(subs$chr[1])
    
    positions <- subs$pos
    
    katLoci = (segInterDist<=kataegis.threshold) # flag specifying if a point is in a peak
    
    if(sum(katLoci)>0) {
      
      start.regions = which(katLoci[-1] & !(katLoci[-(length(katLoci))]) # katLoci breakpoints			
                            | (katLoci[-1] & katLoci[-(length(katLoci))] & segInterDist[-1] != segInterDist[-length(katLoci)] )
      )+1 # endpoints between peaks
      if (katLoci[1]) {start.regions <- c(1, start.regions)}
      
      end.regions = which(!(katLoci[-1]) & katLoci[-(length(katLoci))] #
                          | (katLoci[-1] & katLoci[-(length(katLoci))] & segInterDist[-1] != segInterDist[-length(katLoci)] )
      ) #
      if (katLoci[length(katLoci)]) {end.regions <- c( end.regions, length(katLoci))}
      
      start.regions.init <- start.regions
      end.regions.init <- end.regions 
      
      # handling special cases
      if(length(end.regions)+length(start.regions)>0) {  # if there are any discontinuities in the segmentation at all 						
        if (length(end.regions)==1 & length(start.regions)==0){ 
          start.regions <- 1                                    
        } else if (length(start.regions)==1 & length(end.regions)==0){                                    
          end.regions <- length(positions)                                    
        } else if ((end.regions[1]<start.regions[1])&& (start.regions[length(start.regions)]>end.regions[length(end.regions)])) {
          # starts and ends are the same length, but missing both endpoints
          
          start.regions <- c(1,start.regions)
          end.regions <- c(end.regions,  length(positions))
          
        } else if (end.regions[1]<start.regions[1]){
          # starts will be one shorter
          start.regions <- c(1, start.regions)
          
        } else if (start.regions[length(start.regions)]>end.regions[length(end.regions)]){
          # ends will be one shorter
          
          end.regions <- c(end.regions,  length(positions))
        }
        
        if (length(start.regions)!=length(end.regions)) {
          browser()
        }
        
        
        
        # prepare a data structure that will be later filled up
        kat.regions.all <- data.frame(
          chr=subs$chr[1],
          start.bp=rep(NA,length(start.regions)), # start coordinate [bp]
          end.bp=rep(NA,length(start.regions)), # end coordinate [bp]
          length.bp=rep(NA,length(start.regions)), # length [bp]
          number.bps=rep(NA,length(start.regions)),
          number.bps.clustered=rep(NA,length(start.regions)),
          avgDist.bp=rep(NA,length(start.regions)),
          no.samples=rep(NA,length(start.regions)),
          no.del =rep(NA,length(start.regions)),
          no.dup =rep(NA,length(start.regions)),
          no.inv= rep(NA,length(start.regions)),
          no.trn = rep(NA,length(start.regions)),
          firstBp=start.regions,
          lastBp=end.regions                                    )
        
        kat.regions.all <- hotspotInfo(kat.regions.all, subs, segInterDist)
        
        step.segInterDist.left <- rep(NA, length(segInterDist))
        step.segInterDist.left[2:length(segInterDist)] <- segInterDist[2:length(segInterDist)]- segInterDist[1:(length(segInterDist)-1)]       
        step.segInterDist.right <- rep(NA, length(segInterDist))
        step.segInterDist.right[1:(length(segInterDist)-1)] <- segInterDist[1:(length(segInterDist)-1)]- segInterDist[2:(length(segInterDist))]
        
        kat.regions.all$step.left <-  step.segInterDist.left[start.regions]
        kat.regions.all$step.right <-  step.segInterDist.right[end.regions]
        
        
        # run the filters on the regions of increased frequency
        # make sure there are at least kmin samples
        
        if ((!is.null(kat.regions.all)) && (nrow(kat.regions.all)>0)) {
          kat.regions.all <- subset(kat.regions.all, no.samples>=kmin.samples)
        }
        
        
        # make sure there are at least kmin.filter breakpoints
        if (!is.na(kmin.filter)) {
          kat.regions.all <- subset(kat.regions.all, number.bps>=kmin.filter)
        }
        
        
        
        # make sure the p-value is less than somethng
        if ((!is.null(kat.regions.all)) && (nrow(kat.regions.all)>0)) {
          kat.regions.all <- assignPvalues(kat.regions.all, subs, bp.rate=bp.rate)
          kat.regions.all <- subset(kat.regions.all, pvalue<=pvalue.thresh)
          # only keep the hotspots that exceed the theshold
          kat.regions.all <- subset(kat.regions.all, rate.factor>=rate.factor.thresh)
        }  
        
        
        
        
        
        # merge segments if both were found to be peaks
        if (doMerging) {
          if(nrow(kat.regions.all)>1){
            for(r in 2:nrow(kat.regions.all)){
              if (kat.regions.all$lastBp[r-1] == (kat.regions.all$firstBp[r]-1)) {
                # merge two segments
                kat.regions.all$firstBp[r] <- kat.regions.all$firstBp[r-1]
                kat.regions.all$firstBp[r-1] <- NA
                kat.regions.all$lastBp[r-1] <- NA
                kat.regions.all$avgDist.bp[r] <- NA # this will need to be updated as segments are being merged
              }
            }        
          }
          # remove some of the merged segments
          kat.regions.all <- subset(kat.regions.all, !is.na(firstBp) & !is.na(lastBp))
          
          # update the info on hotspots that might have changed when they were merged
          kat.regions.all <- hotspotInfo( kat.regions.all ,  subs, segInterDist)
          kat.regions.all <- assignPvalues(kat.regions.all, subs, bp.rate=bp.rate)
        } # end merging
        
        
        
        
      } # end if there are discontinuities in the segmentation
    } # if there are any points under the inter-mutation distance threshold
    
    kat.regions.all
    
  }
  
  exactPcf <- function(y, kmin=5, gamma, yest) {
    ## Implementaion of exact PCF by Potts-filtering
    ## x: input array of (log2) copy numbers
    ## kmin: Mininal length of plateaus
    ## gamma: penalty for each discontinuity
    N <- length(y)
    yhat <- rep(0,N);
    if (N < 2*kmin) {
      if (yest) {
        return(list(Lengde = N, sta = 1, mean = mean(y), nIntervals=1, yhat=rep(mean(y),N)))
      } else {
        return(list(Lengde = N, sta = 1, mean = mean(y), nIntervals=1))
      }
    }
    initSum <- sum(y[1:kmin])
    initKvad <- sum(y[1:kmin]^2)
    initAve <- initSum/kmin;
    bestCost <- rep(0,N)
    bestCost[kmin] <- initKvad - initSum*initAve
    bestSplit <- rep(0,N)
    bestAver <- rep(0,N)
    bestAver[kmin] <- initAve
    Sum <- rep(0,N)
    Kvad <- rep(0,N)
    Aver <- rep(0,N)
    Cost <- rep(0,N)
    kminP1=kmin+1
    for (k in (kminP1):(2*kmin-1)) {
      Sum[kminP1:k]<-Sum[kminP1:k]+y[k]
      Aver[kminP1:k] <- Sum[kminP1:k]/((k-kmin):1)
      Kvad[kminP1:k] <- Kvad[kminP1:k]+y[k]^2
      bestAver[k] <- (initSum+Sum[kminP1])/k
      bestCost[k] <- (initKvad+Kvad[kminP1])-k*bestAver[k]^2
    }
    for (n in (2*kmin):N) {
      yn <- y[n]
      yn2 <- yn^2
      Sum[kminP1:n] <- Sum[kminP1:n]+yn
      Aver[kminP1:n] <- Sum[kminP1:n]/((n-kmin):1)
      Kvad[kminP1:n] <- Kvad[kminP1:n]+yn2
      nMkminP1=n-kmin+1
      Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)]+Kvad[kminP1:nMkminP1]-Sum[kminP1:nMkminP1]*Aver[kminP1:nMkminP1]+gamma
      Pos <- which.min(Cost[kminP1:nMkminP1])+kmin
      cost <- Cost[Pos]
      aver <- Aver[Pos]
      totAver <- (Sum[kminP1]+initSum)/n
      totCost <- (Kvad[kminP1]+initKvad) - n*totAver*totAver
      
      if (length(totCost)==0 || length(cost)==0) {
        browser()
      }
      if (totCost < cost) {
        Pos <- 1
        cost <- totCost
        aver <- totAver
      }
      bestCost[n] <- cost
      bestAver[n] <- aver
      bestSplit[n] <- Pos-1
    }
    n <- N
    antInt <- 0
    if(yest){
      while (n > 0) {
        yhat[(bestSplit[n]+1):n] <- bestAver[n]
        n <- bestSplit[n]
        antInt <- antInt+1
      }
    } else {
      while (n > 0) {
        n <- bestSplit[n]
        antInt <- antInt+1
      }
    }
    n <- N  #nProbes 
    lengde <- rep(0,antInt)
    start <- rep(0,antInt)
    verdi <- rep(0,antInt)
    oldSplit  <- n
    antall <- antInt
    while (n > 0) {
      start[antall] <- bestSplit[n]+1
      lengde[antall] <- oldSplit-bestSplit[n]
      verdi[antall] <- bestAver[n]
      n <- bestSplit[n]
      oldSplit <- n
      antall <- antall-1
    }
    if (yest) {
      return(list(Lengde = lengde, sta = start, mean = verdi, nIntervals=antInt, yhat=yhat))
    } else {
      return(list(Lengde = lengde, sta = start, mean = verdi, nIntervals=antInt))
    }
  }
  
  
  medianFilter <- function(x,k){
    #https://rdrr.io/bioc/copynumber/src/R/medianFilter.r
    n <- length(x)
    filtWidth <- 2*k + 1
    
    #Make sure filtWidth does not exceed n
    if(filtWidth > n){
      if(n==0){
        filtWidth <- 1
      }else if(n%%2 == 0){
        #runmed requires filtWidth to be odd, ensure this:
        filtWidth <- n - 1
      }else{
        filtWidth <- n
      }
    }
    
    runMedian <- runmed(x,k=filtWidth,endrule="median")
    
    return(runMedian)
    
  }
  
  getMad <- function(x,k=25){
    
    #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
    x <- x[x!=0]
    
    #Calculate runMedian
    runMedian <- medianFilter(x,k)
    
    dif <- x-runMedian
    SD <- mad(dif)
    
    return(SD)
  }
  
  calcIntermutDist <- function (subs.type, first.chrom.na = FALSE) {
    
    subs.type.processed <- data.frame()
    for (c in unique(subs.type$chr)) {
      # choose subs from only one chromosome at a time
      
      subs.type.chrom <- subset(subs.type, subset=subs.type$chr==c)
      # sort the subs by position
      subs.type.chrom <- subs.type.chrom [order(subs.type.chrom$position),]
      
      if (first.chrom.na) {
        subs.type.chrom$prevPos <- c(NA,subs.type.chrom$position[1:nrow(subs.type.chrom)-1])
      } else {
        subs.type.chrom$prevPos <- c(0,subs.type.chrom$position[1:nrow(subs.type.chrom)-1])        
      }
      subs.type.chrom$distPrev  <- subs.type.chrom$position -  subs.type.chrom$prevPos
      
      subs.type.processed <- rbind(subs.type.processed,subs.type.chrom)
    }
    
    subs.type.processed$distPrev[subs.type.processed$distPrev==0] <- 1
    subs.type.processed 
  }
  
  
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
    
    input_matrix <- matrix(NA, nrow = 1, ncol = 6,
                            dimnames = list(
                             sample_name,
                             c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
                            )
                           )
    
    if(is.null(indels_class)){
      input_matrix[sample_name,"del.mh.prop"] <- 0
    }else{
      indel_ls  <- summarize_indels(indels_class)
      input_matrix[sample_name,"del.mh.prop"] <- indel_ls$sim_proportions["del.mh.prop"]
    }
    
    if(is.null(snv_df)){
      input_matrix[sample_name,"SNV3"] <- 0
      input_matrix[sample_name,"SNV8"] <- 0
    }else{
      SNV_ls    <- summarize_SNVs_dcSigs(snv_df, SNV_sigs, sample_name, subsample_proportion=0.8)
      input_matrix[sample_name,"SNV3"] <- SNV_ls$SBS3
      input_matrix[sample_name,"SNV8"] <- SNV_ls$SBS8
    }
    
    if(is.null(SV_catalogue)){
      input_matrix[sample_name,"SV3"]  <- 0
      input_matrix[sample_name,"SV5"]  <- 0
    }else{
      SV_ls     <- summarize_SVs(SV_catalogue, SV_sigs, sample_name, nboot = 200)
      input_matrix[sample_name,"SV3"]  <- SV_ls$SV_sigtools_exposures$RefSigR3
      input_matrix[sample_name,"SV5"]  <- SV_ls$SV_sigtools_exposures$RefSigR5
    }
    
    LOH_ls    <- summarize_LOH(seg.data)
    input_matrix[sample_name,"hrd"]  <- unlist(LOH_ls[1])
    
    suppressMessages(
      HRDetect_res <- HRDetect_pipeline(data_matrix = input_matrix, bootstrapHRDetectScores = FALSE)
    )

    boot_results <- c(HRDetect_res$data_matrix, HRDetect_res$hrdetect_output)
    hrdetect_boots <- rbind.data.frame(hrdetect_boots,boot_results )
    names(hrdetect_boots) <- c(colnames(HRDetect_res$data_matrix), paste0(colnames(HRDetect_res$hrdetect_output),".w"))
    
  }   
  cat('\ndone bootstrapping\n')
  
  quantiles <- lapply(hrdetect_boots, function (x) quantile(x, probs = c(5, 50, 95)/100))
  
  return(quantiles)  
}

reformat_signatures <- function(SNV_sigs){
  rownames(SNV_sigs) <- SNV_sigs[,1]
  SNV_sigs <- SNV_sigs[,-1]
  SNV_sigs_t <- as.data.frame(t(SNV_sigs))
  return(SNV_sigs_t)
}

sigTools_formatter <- function(input,sampleName){
  
  cat_list <- list()
  names(input)[1] <- "catalogue"
  colnames(input$catalogue) <- sampleName
  cat_list[[1]] <- input$catalogue
  catalogues <- do.call(cbind,cat_list)
  return(catalogues)
  
}

simpleEventType <- function(gr) {
  #TODO: counts identified events twice
  
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

summarize_indels <- function(indels_class, seed=NULL){
  
  ## Using classified in/del mutations, shuffle counts according to proportions
  
  set.seed(seed)
  
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

summarize_SNVs_dcSigs <- function(snv_df, SNV_sigs, sample_name, subsample_proportion, seed=NULL ){
  
  set.seed(seed)
  
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

summarize_SVs <- function(SV_catalogue, SV_sigs, sampleName , nboot = 200, giniThresholdScaling=2, threshold_percent=1, seed=NULL){
    
  set.seed(seed)
  

  SV_catalogue_reformat <-  sigTools_formatter(SV_catalogue, sampleName = sample_name)
  
  #fit SV signature
  SV_fit_sigs_all <- Fit(catalogues = SV_catalogue_reformat,
                         signatures = SV_sigs,
                         useBootstrap = TRUE,
                         nboot = nboot,
                         exposureFilterType = "giniScaledThreshold",
                         giniThresholdScaling = giniThresholdScaling, 
                         threshold_percent = threshold_percent,
                         randomSeed = seed
                        )
  
  #extract signature exposures
  SV_exposure <- as.data.frame(t(SV_fit_sigs_all$exposures[1,]))
  
  #make SV list
  SV_ls <- list(SV_catalogue$rearr_catalogue,SV_exposure)
  
  names(SV_ls) <- c("SV_sigtools_catalog","SV_sigtools_exposures")
    
  return(SV_ls)
  
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


