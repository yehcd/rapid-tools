#' # LIBRARIES
#' requires commonFunctions.Rmd
#' strawR for HiC: https://github.com/aidenlab/straw/tree/master/R
## ----------------------------------------------------------------------------------
library(strawr)

#' 
#' # FUNCTION - ALL-IN-ONE LOADER
#' Wrapper for functions prepareHicFor4C() & prepareRapidFor4C() intended to aggregate RaPID-seq, 4C/HiC data, and DSB site information to generate a single GRanges object containing all of the necessary information for further processing.
#' 
#' This could be sped up if the Hic/4C loading were simply reduced to the chromosome of interest rather than loading the whole genome.
## ----------------------------------------------------------------------------------
aggregateRapidHic <-
  function(pathToHiC,
           pathToRapid,
           dsbSite,
           resolution = 25000,
           regionSize = 3E6,
           dsbMaskSize = 10 * 15000,
           bgLift = FALSE) {
    # load 4C data slice from .hic datafile source
    grl.hic.4c <-
      prepareHicFor4C(
        pathToHic = pathToHiC,
        dsbSite = dsbSite,
        regionSize = regionSize,
        resolution = resolution,
        dsbMaskSize = dsbMaskSize
      )
    
    # load RaPID-seq data & rescale per 4C 'expected' contact model.
    gr.rapid <- loadBigwigToGranges(pathToRapid)
    
    grl.rapid.4c <-
      prepareRapidFor4C(
        gr.rapid = gr.rapid,
        gr.4c.expected = grl.hic.4c$"expected",
        dsbSite = dsbSite,
        regionSize = regionSize,
        resolution = resolution,
        dsbMaskSize = dsbMaskSize,
        bgLift = bgLift
      )
    
    # aggregate data from loading into a single GRanges object
    gr.aggregated <- grl.hic.4c$expected %>% granges()
    gr.aggregated$masked <- grl.hic.4c$expected$masked
    gr.aggregated$midpoint <- grl.hic.4c$expected$midpoint
    gr.aggregated$distance <- gr.aggregated$midpoint %>% abs()
    
    gr.aggregated$rapid.obs <- grl.rapid.4c$observed$score
    gr.aggregated$hic.obs <-  grl.hic.4c$observed$score
    gr.aggregated$expected <- grl.hic.4c$expected$score
    
    return(gr.aggregated)
    
  }


#' 
#' 
#' 
#' # FUNCTION - MAIN
#' 
#' This function handles all-in-one loading of HIC data, specificially:
#'   1. loading of HIC data from .hic datafile
#'   2. slice out 4C view from DSB site
#'   3. rescaling such that maximum value outside of dsbMask region is set to 1
#'   
#' Returns a GRangesList containing data for rescaled 4C EXPECTED & OBSERVED counts of the dsbSite view.
#' 
#' pathToHic         As string. path to .hic datafile. i.e. "ENCFF621AIY.hic"
#' 
#' dsbSite           As string. given as bp position immediately following the DSB site (e.g., for
#'                   gCY201/HBB/J10), the input should be given as "chr11:5,226,985".
#'                   
#' regionSize        BP of region to extract. Given as half-width; default is 3E6 for a total window size 
#'                   of 6Mbps. See makeDsbGranges() for specific details.
#' 
#' resolution        BP resolution to extract from HiC data file. This should match the RaPID-seq binSize! By 
#'                   default this is 25kbps. Valid values given by "strawr:readHicBpResolutions(pathToHiC)"
#'                   
#' dsbMaskSize       Region around DSB to exclude when performing rescaling. Given as half-width; default is                        150kbps for a total window size of 300kbps. See makeDsbGranges() for specific details.
#' 
## ----------------------------------------------------------------------------------
prepareHicFor4C <-
  function(pathToHic,
           dsbSite,
           regionSize = 3E6,
           resolution = 25000,
           dsbMaskSize = 10 * 15000) {
    
    # setup region of interest and dsbMask
    gr.dsbRegion <-
      makeDsbGranges(dsbSite = dsbSite,
                     windowSize = regionSize,
                     resolution = resolution)
    
    gr.dsbMask <-
      makeDsbGranges(dsbSite = dsbSite,
                     windowSize = dsbMaskSize,
                     resolution = resolution)
    
    # extract whole chromosome 4C views for DSB site
    gr.4c.expected <-
      retrieveHicWholeGenome(
        pathDataHic = pathToHiC,
        dsbSite = dsbSite,
        resolution = resolution,
        typeHicMatrix = "expected"
      )
    
    gr.4c.observed <-
      retrieveHicWholeGenome(
        pathDataHic = pathToHiC,
        dsbSite = dsbSite,
        resolution = resolution,
        typeHicMatrix = "observed"
      )
    
    # subset 4C data to region of interest (i.e., regionSize around DSB)
    gr.4c.expected.roi <-
      subset4cToDsb(gr.4c = gr.4c.expected,
                    gr.dsb = gr.dsbRegion,
                    resolution = resolution)
    
    gr.4c.observed.roi <-
      subset4cToDsb(gr.4c = gr.4c.observed,
                    gr.dsb = gr.dsbRegion,
                    resolution = resolution)
    
    
    # rescale 4C data locally by excluding dsbMaskSize region
    grl.4c.scaled <-
      rescaleByDsbMask2(
        gr.4c.roi.exp = gr.4c.expected.roi,
        gr.4c.roi.obs = gr.4c.observed.roi,
        gr.dsbMask = gr.dsbMask
      )
    
    return(grl.4c.scaled)
  }

#' 
#' 
#' This is a wrapper function that takes RaPID-seq data and rescales it for comparison to 4C data by estimating the RaPID-seq "pure proximity" values values on 4C "EXPECTED" values. Requires results from prepareHicFor4C(). 
#' 
#' Simply calls estimateRapid4cBackground() & rescaleByDsbMask2().
#' 
#' gr.rapid          output of loadBigwigToGranges()
#' 
#' gr.4c.expected    output of prepareHicFor4C (GRangesList); specifically GRanges dataset named "expected". 
#' 
#' dsbSite           As string. given as bp position immediately following the DSB site (e.g., for
#'                   gCY201/HBB/J10), the input should be given as "chr11:5,226,985".
#'                   
#' regionSize        BP of region to extract. Given as half-width; default is 3E6 for a total window size 
#'                   of 6Mbps. See makeDsbGranges() for specific details.
#' 
#' resolution        BP resolution to extract from HiC data file. This should match the RaPID-seq binSize!
#'                   By default this is 25kbps. Valid values given by
#'                   "strawr:readHicBpResolutions(pathToHiC)"
#'                   
#' dsbMaskSize       Region around DSB to exclude when performing rescaling. Given as half-width; default is                        150kbps for a total window size of 300kbps. See makeDsbGranges() for specific details.
#' 
## ----------------------------------------------------------------------------------
prepareRapidFor4C <-
  function(gr.rapid,
           gr.4c.expected,
           dsbSite,
           regionSize = 3E6,
           resolution = 25000,
           dsbMaskSize = 10 * 15000,
           bgLift = FALSE) {
    
    
    # setup region of interest and dsbMask
    gr.dsbRegion <-
      makeDsbGranges(dsbSite = dsbSite,
                     windowSize = regionSize,
                     resolution = resolution)
    
    gr.dsbMask <-
      makeDsbGranges(dsbSite = dsbSite,
                     windowSize = dsbMaskSize,
                     resolution = resolution)
    
    # subset RaPID-seq data to region of interest around DSB site
    gr.rapid.roi <-
      subset4cToDsb(gr.4c = gr.rapid,
                    gr.dsb = gr.dsbRegion,
                    resolution = resolution)
    
    
    #### EXPERIMENTAL
    # TODO: testing value lift to eliminate negative values
    
    # if negative, then lift floor to "zero"; 
    
    # but if positive; "drop" to zero OR lease as-is?
    # dropping might be useful to normalize high local background.
    # the "danger" of this assumption is if the local 4C.observed is high
    # in the entire rapid.roi due to bona fide chromatin context.
    
    # at sufficently high regionSize (5E6?) this should be unlikely
    
    # perhaps "safer" to lift by globalMedian (whole genome)?
    
    # alternatively, can simply drop negative values without "lift"
    # the assertion being, if negative values (noDSB enriched) are very high
    # then most likely the RaPID-seq signal at the ROI is low (e.g., few DSBs
    
    if (bgLift) {
      rapidMinScore <- min(gr.rapid.roi$score)
      gr.rapid.roi$score <- gr.rapid.roi$score - rapidMinScore
    }

    #### END BLOCK
    
    
    # estimate RaPID-seq background based on 4C expected model for DSB-site 
    gr.rapid.4c.background <- estimateRapid4cBackground(
      gr.rapid.roi = gr.rapid.roi,
      gr.4c.expected = gr.4c.expected,
      gr.dsbMask = gr.dsbMask,
      # gr.dsbRegion = gr.dsbRegion,
      dsbSite = dsbSite,
      dsbScalingRegionSize = 2E6,
      resolution = 25000
    )
    
    # rescale 4C data locally by excluding dsbMaskSize region
    grl.rapid.4c.scaled <-
      rescaleByDsbMask2(
        gr.4c.roi.exp = gr.rapid.4c.background,
        gr.4c.roi.obs = gr.rapid.roi,
        gr.dsbMask = gr.dsbMask
      )
    
    return(grl.rapid.4c.scaled)
  }

#' 
#' # FUNCTIONS - extra
#' 
#' copied from analyzeRaPID(). DO NOT CHANGE. Use to 'fix' output of strawR to fill in missing windows with 0 score.
## ----------------------------------------------------------------------------------
# using GRanges:binnedAverage() to aggregate RaPID-seq FKPM normalized data into
# regularly-sized windows of target width

binRaPID <-
  function(genome,
           tileWidth = 25000,
           gr.RaPID,
           weight = "score") {
    # load genome GRanges and create adjacent tiles at HiC data tileWidth
    gr.genome <-
      GRangesForBSGenome(genome = genome) %>% keepStandardChromosomes(pruning.mode = "coarse")
    
    # make windows at set resolution
    gr.genome.windows <-
      slidingWindows(x = gr.genome, width = tileWidth, step = tileWidth) %>% unlist()
    
    
    # name of quantitative metric is "score" by default
    rlel.RaPID.cov <- coverage(gr.RaPID, weight = weight)
    gr.RaPID.binnedAvg <-
      binnedAverage(bins = gr.genome.windows,
                    numvar = rlel.RaPID.cov,
                    varname = weight)
    
    return(gr.RaPID.binnedAvg)
  }

#' 
#' 
#' Using 4C 'view' from HiC data for the "OBSERVED" dataset, calculate RaPID-values for a ROI to estimate background model for pure-proximity. 
#' 
#' In the case of RaPID-seq, a DSB-region mask should be applied due to known signal depletion due to resection (e.g., DSB+- (30-100)kbps). All inputs given as GRanges objects.
#' 
#' THE RESOLUTION MUST BE THE SAME BETWEEN: (gr.rapid) & (gr.4c.expected); typically, 25kbps.
#' 
#' gr.rapid              output of loadBigwigToGranges() on analyzeRaPID() bigWig outputs at target          
#'                       resolution. 
#' 
#' gr.4c.expected        output of retrieveHicWholeGenome() for a DSB site for typeHicMatrix = "expected"; a
#'                       GRanges representing the 4C 'view' from the DSB-containing bin.
#' 
#' gr.dsbMask            output of makeDsbGranges(); a GRanges that is centered at the DSB and defines the
#'                       region to mask from further analysis.
#' 
#' <!-- gr.dsbRegion          output of makeDsbGranges(); a GRanges representing the region around the DSB to use -->
#' <!--                       in the estimation. ONLY values in this region will be considered. Typically -->
#' <!--                       DSB+-3Mbps. -->
#'                       
#' dsbSite           As string. given as bp position immediately following the DSB site (e.g., for
#'                   gCY201/HBB/J10), the input should be given as "chr11:5,226,985".
#'                   
#' regionSize        BP of region to extract. Given as half-width; default is 3E6 for a total window size 
#'                   of 6Mbps. See makeDsbGranges() for specific details.
#' 
#' resolution        BP resolution to extract from HiC data file. This should match the RaPID-seq binSize!
#'                   By default this is 25kbps. Valid values given by
#'                   "strawr:readHicBpResolutions(pathToHiC)"
#' 
#' Returns all tiles present in gr.rapid, with those not being used in background estimation having the mask value set as TRUE.
## ----------------------------------------------------------------------------------
estimateRapid4cBackground <-
  function(gr.rapid.roi,
           gr.4c.expected,
           gr.dsbMask,
           # gr.dsbRegion,
           dsbSite,
           dsbScalingRegionSize = 2E6,
           resolution = 25000) {
    
    # setup region of interest to consider for rescaling (e.g., DSB +/- 2Mbps)
    gr.dsbScalingRegion <-
      makeDsbGranges(dsbSite = dsbSite,
                     windowSize = dsbScalingRegionSize,
                     resolution = resolution)
    
    
    
    # subset GRanges for 4C.expected data to same as gr.rapid.roi (i.e., gr.dsb)
    gr.4c.expected.roi <-
      getDsbOverlaps(gr.data = gr.4c.expected,
                     gr.dsb = gr.rapid.roi,
                     type = "include")
    
    
    ## calculate scaling factors
    # subset 4C data to dsbScalingRegion
    gr.4c.expected.roi.masked <-
      getDsbOverlaps(gr.data = gr.4c.expected,
                     gr.dsb = gr.dsbScalingRegion,
                     type = "include")
    
    # mask region immediately around DSB
    gr.4c.expected.roi.masked <-
      getDsbOverlaps(gr.data = gr.4c.expected.roi.masked,
                     gr.dsb = gr.dsbMask,
                     type = "exclude")

        
    # subset RaPID-0seq data to dsbScalingRegion
    gr.rapid.roi.masked <-
      getDsbOverlaps(gr.data = gr.rapid.roi,
                     gr.dsb = gr.dsbScalingRegion,
                     type = "include")    
    
    # mask region immediately around DSB
    gr.rapid.roi.masked <-
      getDsbOverlaps(gr.data = gr.rapid.roi.masked,
                     gr.dsb = gr.dsbMask,
                     type = "exclude")
    
    
    # get total score for regions (dsbScalingRegion - dsbMask)
    maskedSum.4c.obs <- sum(gr.4c.expected.roi.masked$score)
    maskedSum.rapid <- sum(gr.rapid.roi.masked$score)
   
    
    # calculate per-tile scaling factor
    gr.4c.expected.roi$scale <-
      (gr.4c.expected.roi$score) / maskedSum.4c.obs
    
    
    # make 4c proximity background based on RaPID-seq counts
    gr.output <- gr.rapid.roi %>% granges()
    gr.output$midpoint <- gr.rapid.roi$midpoint
    gr.output$score <- maskedSum.rapid * gr.4c.expected.roi$scale

      
    # save mask information
    maskHits <-
      findOverlaps(query = gr.dsbMask, subject = gr.output)
    
    gr.output$masked <- FALSE
    gr.output[subjectHits(maskHits)]$masked <- TRUE
    
    return(gr.output)
  }

#' 
#' 
#' 
#' 
#' 
#' # FUNCTIONS - LOAD DATA - HiC
#' 
## ----------------------------------------------------------------------------------
# no apparent way to retrieve genome assembly information via StrawR? Need to manually check header! For Aiden K562 2023 (ENCFF621AIY.hic), this is hg38.
# returns output a GRanges

# TODO: VERIFY that different chr 4C interactions being extracted is actually fully interaction set & not partial!
retrieveHiC <-
  function(subjectChr,
           pathDataHic,
           dsbSite,
           resolution,
           typeHicMatrix,
           genome) {
    
# parse dsbSite to chromosome and bp position mixed type list & remove "chr" for strawR compatibility
    dsb <- parseDsbSiteString(dsbSite)
    
    
    # resolve exact dsb site to tile at given resolution
    # TODO: HiC dataformat appears to be 0-based?
    # TODO: HiC values given as START of bins in units "BP" 
    resolved.dsb <-
      (as.integer(dsb[["position"]] / resolution)) * resolution

    
    # load HiC data via strawR
    # default: (norm = "SCALE"; e.g., "Balanced++")
    hic.k562.chromosome <- straw(
      norm = "SCALE",
      fname = pathDataHic,
      chr1loc = dsb[["chr"]],
      chr2loc = subjectChr,
      unit = "BP",
      binsize = resolution,
      matrix = typeHicMatrix
    )
    
    
    # simulate 4C/5C output from HiC data; basically ignore existence of missing values from sparse matrix
    # extract datapoints from strawR output
    
    
    # TODO: find a 'all-in-one' method to build this correctly
    
    if (dsb[["chr"]] == subjectChr)
    {
      df.hic.x <-
        hic.k562.chromosome[hic.k562.chromosome$x == resolved.dsb, ]
      df.hic.y <-
        hic.k562.chromosome[hic.k562.chromosome$y == resolved.dsb, ]

      # flip xy for 'df.hic.y'
      colnames(df.hic.y) <- c("y", "x", "counts")


      # merge halves, remove duplicated value
      df.hic.xy.merged <-
        rbind(df.hic.x, df.hic.y) %>% arrange(y, x)
      df.hic.xy.merged <-
        df.hic.xy.merged[!duplicated(df.hic.xy.merged),]
      

    } else{
      df.hic.xy.merged <-
        hic.k562.chromosome[hic.k562.chromosome$x == resolved.dsb,]
      
      
      # if no HiC data (rows), format & return empty GRanges
      if ((dim(df.hic.xy.merged)[1] == 0)) {
        gr.hic.roi <-
          GRanges(seqinfo = SeqinfoForBSGenome(genome = genome)) %>% keepStandardChromosomes(pruning.mode = "coarse")
        gr.hic.roi$score <- numeric(0)

        return(gr.hic.roi)
      }
    }
    
    # friendly naming & ordering
    colnames(df.hic.xy.merged) <-
      c("querySite", "subjectSite", "counts")
    df.hic.xy.merged$dsbChr <- dsb[["chr"]]
    df.hic.xy.merged$subjectChr <- subjectChr
    
    df.hic.xy.merged <-
      df.hic.xy.merged[, c("dsbChr", "querySite", "subjectChr", "subjectSite", "counts")]
    

    # prepare df for making GRanges output. 
    # set IRanges starts, rename 'counts' to 'score'
    # TODO: verify correct conversion of HiC positions to GRanges base positions
    df.output <-
      df.hic.xy.merged %>% select(subjectChr, subjectSite, counts)
    df.output$subjectChr <- df.output$subjectChr
    df.output$subjectEnd <- df.output$subjectSite + resolution
    names(df.output)[names(df.output) == "counts"] <- "score"
    

    # for unknown reasons (ENCFF621AIY.hic) has out-of-bounds tiles at the end
    # of chr10; trimmed and suppressed.
    suppressWarnings(
      gr.hic.roi <-
        makeGRangesFromDataFrame(
          df = df.output,
          keep.extra.columns = TRUE,
          seqnames.field = "subjectChr",
          start.field = "subjectSite",
          end.field = "subjectEnd",
          starts.in.df.are.0based = TRUE,
          seqinfo = SeqinfoForBSGenome(genome = genome)
        ) %>% keepStandardChromosomes(pruning.mode = "coarse") %>% trim()
    )
    
    return(gr.hic.roi)
  }


# wrapper for retrieveHiC(); retreives a <4C SLICE> based on DSB site 

# builds GRanges for whole genome (human) normal chromosomes
# does not retrieve ChrY since using K562 data

# only works for hg38 (HiC 2022/2023); will not work correctly with hg19 (HiC 2014)
retrieveHicWholeGenome <- function(pathDataHic,
                                   dsbSite,
                                   resolution,
                                   typeHicMatrix) {
  # hardcoded var
  
  genome = "hg38"
  
# built in 'strawr:pathDataHic()' will retrieve chromosome names that do not have corresponding data (i.e., 'All' and 'Y')
  chr.order <- c(
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM"
  )
  
  gr.output <- lapply(
    chr.order,
    retrieveHiC,
    pathDataHic = pathDataHic,
    dsbSite = dsbSite,
    resolution = resolution,
    typeHicMatrix = typeHicMatrix,
    genome = genome
  ) %>% GRangesList() %>% unlist()
  
  
  # hacky fix to set missing tiles from strawR dump to score value of zero
  # also would re-index tiles if starts & ends were out-of-sync for some reason (THIS SHOULD NOT HAPPEN)
  gr.output<- binRaPID(
    genome = genome,
    tileWidth = resolution,
    gr.RaPID = gr.output,
    weight = "score"
  )
  
  
  return(gr.output)
  
}


#' 
#' 
#' extract GRanges for the region around a dsbSite for HiC data, including re-centering such that x=0 is the set to the DSB tile midpoint. 
#' 
## ----------------------------------------------------------------------------------
subset4cToDsb <- function (gr.4c, gr.dsb, resolution) {
  # subset whole genome 4C data to around region of interest (e.g., DSB site)
  gr.4c.roi <- getDsbOverlaps(gr.data = gr.4c,
                              gr.dsb = gr.dsb,
                              type = "include")
  gr.4c.roi$midpoint <- mid(gr.4c.roi)
  
  # get value for centering DSB tile midpoint to x=0
  dsbTileMidpoint <-
    (as.integer(mid(gr.dsb) / resolution) + 0.5) * resolution
  
  # adjust midpoints
  gr.4c.roi$midpoint <- gr.4c.roi$midpoint - dsbTileMidpoint
  
  return(gr.4c.roi)
}


## ----------------------------------------------------------------------------------
# for 4C dataset rescaling of (OBS & EXP) by (EXP) values
rescaleByDsbMask2 <- function(gr.4c.roi.exp, gr.4c.roi.obs, gr.dsbMask) {
  # mask DSB region from ROI
  gr.4c.roi.exp.masked <- getDsbOverlaps(gr.data = gr.4c.roi.exp,
                                     gr.dsb = gr.dsbMask,
                                     type = "exclude")
  
  
  # get scaling factor & rescale to max = 1 (excluding masked region)
  maxScore <- gr.4c.roi.exp.masked$score %>% max()
  
  gr.4c.roi.exp$score <- gr.4c.roi.exp$score / maxScore
  
  # use same scaling factor for "observed" data
  gr.4c.roi.obs$score <- gr.4c.roi.obs$score / maxScore
  
  
  # save mask information
  maskHits <-
    findOverlaps(query = gr.dsbMask, subject = gr.4c.roi.exp)
  
  gr.4c.roi.exp$masked <- FALSE
  gr.4c.roi.exp[subjectHits(maskHits)]$masked <- TRUE
  
  gr.4c.roi.obs$masked <- FALSE
  gr.4c.roi.obs[subjectHits(maskHits)]$masked <- TRUE  
  
  # aggregate and return
  grl.output <- GRangesList(gr.4c.roi.exp, gr.4c.roi.obs)
  names(grl.output) <- c("expected", "observed")
  return(grl.output)
}


#' 
#' 
#' 
#' 
#' # FUNCTIONS - PLOTTING
#' 
#' customized plotting for HiC data only
## ----------------------------------------------------------------------------------
plotGenomicRegionHic <-
  function(xRange = c(-2E6, 2E6),
           yRange = c(-1, 1),
           yBreaks = c(-1, -0.5, 0, 0.5, 1),
           xBreaks = c(-2E6, -1E6, 0 , 1E6, 2E6),
           df.data,
           outputName) {
    plot <- ggplot(data = df.data,
                   aes(x = midpoint,
                       y = score)) + geom_hline(yintercept = 0,
                                                color = "gray40",
                                                linetype = "dotted") + geom_vline(xintercept = 0,
                                                                                  color = "gray40",
                                                                                  linetype = "dotted") + geom_line(aes(x = midpoint, y = score), color = "#D81B60") +
      theme(
        axis.text.y = element_text(margin = margin(r = 10)),
        axis.text.x = element_text(margin = margin(t = 10))
      )  + coord_cartesian(
        xlim = xRange,
        ylim = yRange,
        expand = TRUE
      ) +  scale_x_continuous(breaks = xBreaks)  + scale_y_continuous(breaks = yBreaks) + xlab("distance from DSB (bps)") + ylab("norm. RaPID score") + theme_classic()
    
    
    return(plot)

  }

#' 
#' generic basic plotting function for HiC & RaPID data
## ----------------------------------------------------------------------------------
plotGenomicRegionDual <-
  function(xRange = c(-2E6, 2E6),
           yRange = c(-1, 1),
           yBreaks = c(-1,-0.5, 0, 0.5, 1),
           xBreaks = c(-2E6,-1E6, 0 , 1E6, 2E6),
           df.upper,
           df.lower) {
    plot <- ggplot(data = df.upper,
                   aes(x = midpoint,
                       y = score)) + geom_area(aes(x = midpoint, y = score), fill = "#D81B60") + geom_area(aes(x = df.lower$midpoint, y = -df.lower$score), fill = "#1E88E5") + #, alpha = 0.4) +
      theme(
        axis.text.y = element_text(margin = margin(r = 10)),
        axis.text.x = element_text(margin = margin(t = 10))
      )  + coord_cartesian(xlim = xRange,
                           ylim = yRange,
                           expand = TRUE) +  scale_x_continuous(breaks = xBreaks)  + scale_y_continuous(breaks = yBreaks) + xlab("distance from DSB (bps)") + ylab("scaled RaPID score") + geom_hline(yintercept = 0,
                                                                                                                                                                                                     color = "gray20") + geom_vline(xintercept = 0,
                                                                                                                                                                                                                                       color = "gray20",
                                                                                                                                                                                                                                       linetype = "dotted") + theme_classic()
    
    
    return(plot)

  }


#' 
#' 
#' 
