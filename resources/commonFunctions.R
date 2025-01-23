#' # OVERVIEW
#' This notebook is for consolidating functions used for further RaPID-seq plotting and analysis into a single place to ensure consistency of data processing
#' 
#' 13 March 2024 - Created
#' 
#' 
#' # LIBRARIES
## -------------------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges) 
library(matrixStats)

#' 
#' # FUNCTIONS - GENERIC USE
#' 
#' Splits a given dsbSite (example input: chr11:5,226,985) into list components chr and position.
## -------------------------------------------------------------------------------------
parseDsbSiteString <- function(dsbSite){
  dsbSite <- gsub(pattern = ",", replacement = "", x = dsbSite)
  
  dsbChr <-
    strsplit(dsbSite, ":") %>% unlist() %>% nth(1) %>% as.character()
  dsbPosition <-
    strsplit(dsbSite, ":") %>% unlist() %>% nth(2) %>% as.integer()
  
  dsbParsed <- list(dsbChr, dsbPosition)
  dsbParsed <- setNames(dsbParsed, c("chr", "position"))
  
  return(dsbParsed)
}

#' 
#' converts DSB location to specific tile at a given resolution. 
#' 
#' dsbParsed       output of parseDsbSiteString()
#' resolution      binSize used in analyzeRaPID()
#' 
#' returns start of GRanges tile at given resolution 
## -------------------------------------------------------------------------------------
resolveDsb <- function(dsbParsed, resolution) {
  
  # simple calculations of positions
  start <-
    (as.integer(dsbParsed[["position"]] / resolution)) * resolution + 1
  
  end <- start + resolution - 1
  midpoint <- start + as.integer(resolution/2) - 1  
  
  # construct output
  output <- c(start, midpoint, end)
  names(output) <- c("start", "midpoint", "end")
    
  return(output)
}


#' 
#' 
#' 
#' Intended for making GRanges for a DSB site in hg38. The DSB position is recentered to the containing
#' tile midpoint. The DsbGranges() region is defined as DSB-containing tile +/- windowSize INCLUSIVE of
#' the DSB-containing tile.
#' 
#' dsbSite       DSB site, given as 1-based position immediately following the DSB 
#'               0-based position. For gCY201 (bp following DSB: chr11:5,226,985).
#' 
#' windowSize    bps in EACH direction to make GRanges window. 
#'               i.e., input 3E6 = 6E6 total size
#'               
#' resolution    bps of each tile in the dataset (e.g., 25kbps)
## -------------------------------------------------------------------------------------
makeDsbGranges <-
  function(dsbSite,
           windowSize = 3E6,
           resolution = 25000) {
    
    # parse DSB site
    dsbParsed <- parseDsbSiteString(dsbSite)
    
    # get tile containing DSB
    resolvedDsb <-
      resolveDsb(dsbParsed = dsbParsed, resolution = resolution)
    resolveWindowSize <- ceiling(windowSize / resolution) - 1
    
    # make GRanges
    gr.dsb <-
      GRanges(
        seqnames = dsbParsed$chr,
        ranges = IRanges(
          start = resolvedDsb["start"] - resolveWindowSize * resolution,
          end = resolvedDsb["end"] + resolveWindowSize * resolution
        )
      )
    
    return(gr.dsb)
  }

#' 
#' 
#' 
#' 
#' # FUNCTIONS - LOAD DATA - RAPID-seq
#' Load bigWig data and apply seqlevels/info for hg38
#' 
#' bigWigFile    path to bigWig file to load (assumes hg38)
## -------------------------------------------------------------------------------------
loadBigwigToGranges <- function(bigWigFile) {
  gr.data <- rtracklayer::import(con = bigWigFile, format = "bigWig")
  
  # force seqlevels consistency
  seqlevels(gr.data, pruning.mode="coarse") <-
    seqlevels(BSgenome.Hsapiens.UCSC.hg38)
  
  seqinfo(gr.data) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
  
  # drop nonstandard chromosomes
  gr.data <- gr.data %>% keepStandardChromosomes(pruning.mode = "coarse")

  
  return(gr.data)
}

#' 
#' 
#' Wrapper function to subset a RaPID-seq GRanges based on a defined ROI (e.g., DSB region)
#' gr.rapid      GRanges from which to extract DSB region from
#' gr.dsb        GRanges from makeDsbGranges()
#' type          what GRanges to return:
#'                   "include" = tiles in gr.rapid overlapping gr.dsb
#'                   "exclude" = tiles in gr.rapid NON-overlapping gr.dsb
#' 
#' resolution    binSize used in data (e.g., 25kbps)
## -------------------------------------------------------------------------------------
getDsbOverlaps <- function(gr.data, gr.dsb, type = "include") {
  # simple check of input
  if (!type %in% c("include", "exclude")) {
    stop("Invalid 'type' value")
  }

  
  # subsetting by gr.dsb & type
  if (type == "include") {
    hits <- findOverlaps(query = gr.dsb, subject = gr.data)
    gr.output <- gr.data[subjectHits(hits)]
  } else if (type == "exclude"){
    hits <- findOverlaps(query = gr.dsb, subject = gr.data)
    gr.output <- gr.data[-subjectHits(hits)]
  }

  return(gr.output)
}


#' 
#' 
#' 
#' Used for rebinning/sizing/centering coverage data from bigWigs via deepTools bamCoverage; deepTools merges adjacent bins, which must be split into fixed-width bins for easy plotting in ggplot2 when using bin midpoints.
#' 
#' gr.chip         output of loadBigwigToGranges(). Intended to originate from deepTools bamCoverage
#' gr.dsb          output of makeDsbGranges(). Granges covering DSB site or other ROI.  
#' tileWidth       bp size of new bins to output
## -------------------------------------------------------------------------------------
getDsbOverlapsBinned <- function(gr.chip, gr.dsb, tileWidth = 250) {
  hits <- findOverlaps(query = gr.dsb, subject = gr.chip)
  gr.output <- gr.chip[subjectHits(hits)]
  
  
  ## remap/bin coverage using GRanges functions
  # make tiles of target size
  gr.tiles <-
      slidingWindows(x = gr.dsb, width = tileWidth, step = tileWidth) %>% unlist()
  
  # name of quantitative metric is "score" by default
  rlel.data.cov <- coverage(gr.chip, weight = "score") 
  
  gr.output <-
    binnedAverage(bins = gr.tiles,
                  numvar = rlel.data.cov,
                  varname = "score")
  
  # get midpoint values
  # convert positions to DSB site as position 0
  gr.output$midpoint <- mid(gr.output) - mid(gr.dsb)

  
  
  return(gr.output)
}

#' 
#' 
#' 
#' # FUNCTIONS - NORMALIZE DATA PER SOME FACTOR
#' 
#' Normalization by rescaling using ONLY maximal global amplitude (max, min); not 0-1 rescaling. This might be considered [(-1)to(+1)] rescaling, where 1 is defined by the single largested value in either direction (typically + in standard RaPID-seq).
#' 
#' Used to make Figure 1 data for showing HBB/gCY201 exemplar data when averaged by this method is highly consistent.
#' 
#' **TODO: change to use with loadBigwigToGranges() to standardize data handling process?
#' 
#' bigWigFile    path to bigWig file written out by analyzeRaPID
## -------------------------------------------------------------------------------------
normalizeRapidGlobally3 <- function(bigWigFile){
  gr.rapid <- rtracklayer::import(con = bigWigFile, format = "bigWig")
  
  # get maximum magnitudes
  max <- gr.rapid$score %>% max() %>% abs()
  min <- gr.rapid$score %>% min() %>% abs()
  
  # get scaling factor
  scale <- max(max, min)
  
  gr.rapid$score.norm <- gr.rapid$score/scale

  return(gr.rapid)
}

#' 
#' 
#' Rescaling a RaPID-seq ROI locally, again using unidirectional scaling as in normalizeRapidGlobally3(). 
#' 
#' gr.roi          output of getDsbOverlaps(). Genome-wide RaPID-seq data subsetted to DSB region.
## -------------------------------------------------------------------------------------
normalizeRapidLocally <- function(gr.roi){

  # get maximum magnitudes
  max <- gr.roi$score %>% max() %>% abs()
  min <- gr.roi$score %>% min() %>% abs()
  
  # get scaling factor
  scale <- max(max, min)
  
  # gr.roi$score.norm <- gr.roi$score/scale
  gr.roi$score<- gr.roi$score/scale
  
  return(gr.roi)
}

#' 
#' 
#' modified from rapidTools (0.0.3)
#' * genome can be UCSC name (e.g., "hg38"), GRanges object, or Seqinfo definition
#' 
## -------------------------------------------------------------------------------------
binRaPID <-
  function(genome,
           tileWidth,
           # gatcTileMinimum = 0,
           gr.RaPID,
           weight = "score") {
    # load genome GRanges and create adjacent tiles at HiC data tileWidth
    # gr.genome <- GRangesForUCSCGenome(genome = genome)
    
    if (is(genome, "character")) {
      gr.genome <-
        GRangesForBSGenome(genome = genome) %>% keepStandardChromosomes(pruning.mode = "coarse")
    } else if (is(genome, "GenomicRanges")) {
      gr.genome <- genome 
    } else if (is(genome, "Seqinfo")) {
      gr.genome <- genome %>% GRanges() %>% GenomicRanges::reduce()
    }

    
    
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
#' 
#' 
#' For calculating total score in region of interest, derived from
#' 'step5_summarizePlasmid.Rmd'
## -------------------------------------------------------------------------------------
calcTotalScore <- function(gr.rapid, gr.roi = NULL, gr.gatcMask = NULL){
  # required inputs
  stopifnot(class(gr.rapid) == "GRanges")
  stopifnot(class(gr.roi) == "GRanges")
  
  # subset gr.rapid to gr.roi (region of interest)
  overlapsRoi <- findOverlaps(query = gr.roi, subject = gr.rapid)
  gr.rapid <- gr.rapid[subjectHits(overlapsRoi)] 
  
  
  # if mask GRanges definition is provided, apply mask
  if (!is.null(gr.gatcMask)) {
    # find masked tiles
    overlaps <- findOverlaps(query = gr.gatcMask, subject = gr.rapid)
    
    # get non-masked tiles
    gr.rapid <- gr.rapid[-subjectHits(overlaps)]
  }
  
  # calculate score-width (by kbps)
  gr.rapid$scoreWidth <- gr.rapid$score * width(gr.rapid) / 1E3
  
  # get total score
  totalScore <- gr.rapid$scoreWidth %>% sum()
  
  return(totalScore)
}


#' # PLOTTING FUNCTIONS
#' 
#' Standardized general plot function with easy axis modification and etc. Data from GRanges (e.g, normalizeRapidLocally()) needs to be converted to dataframe first. Midpoints of GRanges tiles also need to be already calculated.
#' 
#' Input data should be CENTERED such that DSB site = position 0. 
#' 
#' xRange    ggplot2; +/-bps from DSB (position 0) to plot; this should be less than limits of data            in df.data to make use of edge-to-edge plotting with ggplot2/coord_cartesian()
#' 
#' yRange    ggplot2; +/-bps from y-axis 0 to plot
#' 
#' xBreaks   ggplot2; axis tick mark positions 
#' yBreaks   ggplot2; axis tick mark positions 
#' 
#' df.data   output of normalizeRapidLocally() or getDsbOverlaps() that has the following additional           conversions:
#'               midpoint defined as GRanges interval midpoint WITH evenly spaced interval widths
#'               x-axis recentered such that DSB is at position 0 (midpoint)
#'               converted to data.frame
#'               
#' In general, probably should copy this plotting function and adjust as necessary for each individual plotting use-case 
#'     
## -------------------------------------------------------------------------------------
plotGenomicRegion <-
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
                                                                                  linetype = "dotted") + geom_line(aes(x = midpoint, y = score), color = "darkred") +
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
#' 
#' 

