## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)

library(dplyr)
library(tidyr)

library(grid)
library(rtracklayer)
library(GenomicRanges)
library(strawr)

library(png)
library(ggplot2)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# given that DSB is 'in-between' bp positions, provide as bp position immediately after cut site  for consistency.
# output in 'GRanges' style format
# perhaps should convert "dsbSite" immediately to GRanges & feed this directly to other functions?
parseDsbSite <- function(dsbSite) {
  dsbSite <- strsplit(x = dsbSite, split = ":") %>% unlist()
  dsbChr <- dsbSite[1] %>% as.character()
  dsbPosition <-
    gsub(pattern = ",",
         replacement = "",
         x = dsbSite[2]) %>% as.integer()
  
  parsedDsb <- list(dsbChr, dsbPosition)
  names(parsedDsb) <- c("chr", "position")
  
  return(parsedDsb)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#load GATC-intervals GFF file
load.gatc.intervals <- function(intervals.file.path) {
  gatc.intervals <-
    read.table(
      gzfile(intervals.file.path),
      header = FALSE,
      sep = "\t",
      col.names = c(
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute"
      ),
      stringsAsFactors = FALSE
    )
  
  #fix formatting
  gatc.intervals$attribute <-
    sub("tileID=", "", gatc.intervals$attribute)
  
  gatc.intervals$tileID <-
    sub("chr[[:alnum:]]+_", "", gatc.intervals$attribute) %>% as.integer()
  
  gatc.intervals$length.bp <-
    gatc.intervals$end - gatc.intervals$start + 1
  
  gatc.intervals 
}


# short function to merge GFF & HTSeq data & calculate FKPMs
calc.DamID.FKPM <-
  function(data.HTseq, gatc.map) {
    
    # merge & calculate FKPM
    # merging of HTSeq with GFF file will exclude the unmapped HTSeq reads (e.g. "__ambiguous", etc.)
    
    data.merged <- merge(
      x = data.HTseq,
      y = gatc.map,
      by.x = "attribute",
      by.y = "attribute"
    )
    
    
    # FPKM calculation
    data.merged$FKPM <-
      (data.merged$count) / ((sum(data.merged$count) / 1E6) * (data.merged$length.bp /
                                                                 1000))
    
    
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
    
    # formatting output
    data.processed <-
      data.merged %>% select("seqname",
                             "start",
                             "end",
                             "tileID",
                             "length.bp",
                             "count",
                             "FKPM")
    
    # seems to output as lexographic for 'tileID' column without this reordering
    data.processed <-
      data.processed[order(match(data.processed$seqname, chr.order),
                           data.processed$tileID),]
    
    return(data.processed)
  }


#### perform Lambda genomoe spike-in correction and background subtraction
#
# This spike-in correction method asserts that spike-in correction should only be made between a given signal/background
# pair that come from the same Nucleofection batch, and the same library preparation batch
#

calc.DamID.enrichmentPaired <-
  function(gatc.map,
           data.signal,
           data.background,
           lambda.signal,
           lambda.background) {
    
    # calculate lambda factors
    # get lambda reads as % of all reads
    lambda.ratio.signal <-
      (lambda.signal / sum(data.signal$count))
    
    lambda.ratio.background <-
      (lambda.background / sum(data.background$count))
    
    # ratio of %lambda for (background / signal)
    lambda.scaling <- (lambda.ratio.background / lambda.ratio.signal)
    
    
    # calculate FKPM for each GATC.tile
    
    signal.processed <-
      calc.DamID.FKPM(data.signal, gatc.map)
    
    background.processed <-
      calc.DamID.FKPM(data.background, gatc.map)
    
    # TODO: figure out if need to do some other way to get proper mapping
    output <- signal.processed %>% select("seqname",
                                          "start",
                                          "end",
                                          "tileID",
                                          "length.bp", )
    
    
    # calculate log2.FKPM or FKPM.diff as indicated with lambda scaling
    
    output$log2.FKPM <-
      log((signal.processed$FKPM * lambda.scaling + 1) / (background.processed$FKPM + 1),
          base = 2)
    
    
    # calculate FKPM difference
    output$diff.FKPM <-
      (signal.processed$FKPM * lambda.scaling) - (background.processed$FKPM)
    
    
    #copy over counts data to output
    
    output$signal.count <-
      signal.processed$count #/ lambda.ratio.signal
    output$background.count <-
      background.processed$count #/ lambda.ratio.background
    
    #copy FPKM values to output
    output$signal.FKPM <- signal.processed$FKPM
    output$background.FKPM <- background.processed$FKPM
    
    return(output)
  }

#
# using GRanges:binnedAverage() to aggregate RaPID-seq FKPM normalized data into 
# regularly-sized windows of target width 
#
binRaPID <-
  function(genome,
           tileWidth,
           gr.RaPID,
           weight = "score") {
    # load genome GRanges and create adjacent tiles at HiC data tileWidth
    # gr.genome <- GRangesForUCSCGenome(genome = genome)
    gr.genome <-
      GRangesForBSGenome(genome = genome) %>% keepStandardChromosomes(pruning.mode = "coarse")
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


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pathToDefinitionFile format 
#
# col 1   prefix for output files, binSize will be appended to each output
# col 2   filenames for signal dataset (i.e., +DSB)
# col 3   value for Lambda spike-in count for signal
# col 4   filenames for background dataset (i.e., +DSB)
# col 5   value for Lambda spike-in count for background
# col 6   binSizes to use & output datafiles; semicolon-separated string of values
#

#' @export
analyzeRaPID <-
  function(gatcMapFile = NULL,
           definitionFile = NULL,
           htseqDataDirectory = NULL,
           genomeName = NULL,
           outputDirectoryName = "resultsAnalyzeRaPID",
           logging = TRUE) {
    
    # make running log
    if (logging) {
      logfile <-
        paste0("analyzeRaPID_", Sys.time() %>% as.numeric(), ".log.txt")
      
      sink(file = logfile,
           append = FALSE,
           split = TRUE)
    }
    
    # minimal checking of inputs
    stopifnot(
      !is.null(gatcMapFile),
      !is.null(definitionFile),
      !is.null(htseqDataDirectory),
      !is.null(genomeName),
      !is.null(outputDirectoryName)
    )
    
    version <- "0.0.0.3"
    
    # extract sample information pairs
    df.sampleDefinitions <-
      read.csv(
        file = definitionFile,
        header = TRUE,
        sep = "\t",
        colClasses = c(
          'character',
          'character',
          'numeric',
          'character',
          'numeric',
          'character'
        )
      )
    
    # console printout checks
    cat("\n**********")
    cat("\nRaPID-seq analysis\t\t ver:", version)
    cat("\nstarted system date & time:\t", Sys.time() %>% as.character())
    cat("\n\n**********")
    
    cat("\nworking directory:\t", getwd())
    cat("\nGATC map file:\t\t\t", gatcMapFile)
    cat("\nsample definition file:\t\t", definitionFile)
    cat("\nspecified genome:\t", genomeName)
    cat("\n\n**********\n")
    
    # make output directory folder
    cat("\n...making output directory:\t", outputDirectoryName)
    if (!dir.exists(outputDirectoryName)) {
      dir.create(outputDirectoryName)
    } 
    
    # load GATC map GFF file
    cat("\n...loading GATC map file\n")
    df.gatcMap <-
      load.gatc.intervals(gatcMapFile)
    
    
    cat("\n...processing data pairs\n")
    cat("\n**********")
    
    # perform RaPID-seq (signal - background) processing for all dataset pairs at specified resolutions
    for (i in 1:nrow(df.sampleDefinitions)) {
      
      # parse sample pair information
      outputNamePrefix <- df.sampleDefinitions[i, 'outputNamePrefix']
      
      pathToSignal <-
        paste0(htseqDataDirectory, df.sampleDefinitions[i, 'signalDatafileName'])
      
      lambdaCountSignal <-
        df.sampleDefinitions[i, 'signalLambdaCount']
      
      pathToBackground <-
        paste0(htseqDataDirectory, df.sampleDefinitions[i, 'backgroundDatafileName'])
      
      lambdaCountBackground <-
        df.sampleDefinitions[i, 'backgroundLambdaCount']
      
      binSizes <-
        df.sampleDefinitions[i, 'binSizes'] %>% strsplit(";") %>% unlist() %>% as.numeric()
      
      
      # console printout checks
      cat("\n\nstarted processing pair in row:\t", i, "of", nrow(df.sampleDefinitions))
      cat("\nsystem date & time:\t", Sys.time() %>% as.character())
      
      cat("\n\n\tsignal datafile:\t\t", pathToSignal)
      cat("\n\tsignal lambda counts:\t\t", lambdaCountSignal)
      cat("\n\tbackground datafile:\t\t", pathToBackground)
      cat("\n\tbackground lambda counts:\t", lambdaCountBackground)
      
      
      # load datasets
      cat("\n\n\t...loading datafiles")  
      
      df.signal <- read.table(
        gzfile(pathToSignal),
        header = FALSE,
        sep = "\t",
        col.names = c("attribute", "count"),
        stringsAsFactors = FALSE
      )
      
      df.background <- read.table(
        gzfile(pathToBackground),
        header = FALSE,
        sep = "\t",
        col.names = c("attribute", "count"),
        stringsAsFactors = FALSE
      )
      
      
      # perform raw DamID differential enrichment calculation     
      # probably could make memory footprint smaller somehow
      
      cat("\n\n\t...calculating RaPID-seq enrichment")  
      
      # changed from absolute Lambda values to +/-DSB relative pair comparison
      df.DamID.enrichment <- calc.DamID.enrichmentPaired(
        gatc.map = df.gatcMap,
        data.signal = df.signal,
        data.background = df.background,
        lambda.signal = lambdaCountSignal,
        lambda.background = lambdaCountBackground
      )
      
      df.DamID.enrichment <-
        select(df.DamID.enrichment, seqname, start, end, diff.FKPM)
      
      colnames(df.DamID.enrichment) <-
        c("seqname", "start", "end", "score")
      
      
      gr.DamID.enrichment <- makeGRangesFromDataFrame(
        df = df.DamID.enrichment,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "seqname",
        start.field = "start",
        end.field = "end",
        seqinfo = Seqinfo(genome = genomeName)
      ) %>% keepStandardChromosomes(pruning.mode = "coarse")
      
      
      # write bigWig output for unbinned data
      
      outputFile <-
        paste0(outputDirectoryName,
               "/",
               outputNamePrefix, "_raw.bw")
      
      cat(paste0("\n\n\t...writing bigWig output for unbinned data"))    
      rtracklayer::export(object = gr.DamID.enrichment,
                          con = outputFile,
                          format = "bigWig")      
      
      
      # bin raw DamID differential enrichment values to target resolutions via GRanges::binnedAverage
      
      cat("\n\n\tbinning RaPID-seq to resolution:")
      
      for (j in 1:length(binSizes)) {
        
        resolution <- binSizes[j]
        outputFile <-
          paste0(outputDirectoryName,
                 "/",
                 outputNamePrefix,
                 "_binned-" ,
                 resolution/1000,
                 "kbps.bw")
        
        
        cat(paste0("\n\n\t\t...",resolution, " bps"))  
        
        # bin data
        gr.binnedRaPID <- binRaPID(
          genome = genomeName,
          tileWidth = resolution,
          gr.RaPID = gr.DamID.enrichment,
          weight = "score"
        )
        
        cat("\n\t\t...writing bigWig outputs to:\t",
            paste0("./", outputFile))
        
        # write bigWig output
        rtracklayer::export(
          object = gr.binnedRaPID,
          con = outputFile,
          format = "bigWig"
        )
      }
      
      # record completed pair processing
      cat("\n\n\tcompleted processing pair in row:\t", i)
      cat("\n\tsystem date & time:\t", Sys.time() %>% as.character())
      cat("\n")
    }
    
    cat("\n\n**********")
    cat("\ncompleted processing all pairs")
    cat("\nsystem date & time:\t", Sys.time() %>% as.character())
    
    #return()
  }

