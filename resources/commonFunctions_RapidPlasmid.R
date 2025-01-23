
library(ggplot2)

#' 
#' 
#' # STEP 1 - FKPM calculation
#' Common function
## ----------------------------------------------------------------------------------
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
  
  return(gatc.intervals)
}

#' 
#' FOR HG38-ALIGNED HT-SEQ PROCESSING
## ----------------------------------------------------------------------------------
# short function to merge GFF & HTSeq data & calculate FKPMs

# REF for FKPM calculation:
# https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

# aggregate data from DamID for reads aligning to hg38 and reads aligning to plasmids

calc.DamID.FKPM.Hg38.withPlasmids <-
  function(data.HTseq.hg38,
           gatc.map.hg38,
           data.HTseq.plasmid,
           gatc.map.plasmid,
           dsbMask) {
    # set chromosome order for plasmids
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
    
    # merging of HTSeq with GFF file
    # will exclude the unmapped HTSeq reads (e.g. "__ambiguous", etc.)
    
    data.merged.hg38 <- merge(
      x = data.HTseq.hg38,
      y = gatc.map.hg38,
      by.x = "attribute",
      by.y = "attribute"
    )
    
    data.merged.plasmid <- merge(
      x = data.HTseq.plasmid,
      y = gatc.map.plasmid,
      by.x = "attribute",
      by.y = "attribute"
    )
    
    
    # ADDED 9 October 2024
    #
    # perform before calculating FKPM to ensure proper denominator rebalancing
    # only perform is a mask is defined
    
    if (!is.null(dsbMask)) {
      # parse 'dsbMask' to (1) seqname; (2) start; (3) end
      dsbParsed = dsbMask %>% strsplit(split = "[:-]") %>% unlist()
      
      
      # select tiles to mask & remove if ANY part of GATC tile overlaps with mask
      maskIndex <- (data.merged.hg38$seqname == dsbParsed[1]) &
        (data.merged.hg38$end >= as.integer(dsbParsed[2])) &
        (data.merged.hg38$start <= as.integer(dsbParsed[3]))
      
      
      # remove masked GATC tiles from dataset and save
      data.merged.hg38 <- data.merged.hg38[!maskIndex,]
    }
    
    
    
    
    # extract totals values (defined at reads to hg38 + reasd to plasmid)
    totalReads <-
      sum(data.merged.hg38$count) + sum(data.merged.plasmid$count)
    # totalReads.plasmid <- data.merged.plasmid$count %>% sum()
    # return(totalReads)
    
    
    # calculate modified FKPM value as follows:
    #   (fragments per 1kbps mapped to plasmid GATC tiles)
    #   per (1E6 reads mapped to hg38)
    
    data.merged.hg38$FKPM <-
      (data.merged.hg38$count) / ((totalReads / 1E6) * (data.merged.hg38$length.bp / 1000))
    
    
    
    # formatting output
    data.processed <-
      data.merged.hg38 %>% select("seqname",
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



calculateRapidPlasmidFkpmForHg38 <-
  function(pathToSampleInfo,
           genomeName = NULL,
           pathToGatcHg38,
           pathToGatcPlasmids,
           inputDataDir,
           outputDir,
           dskMask = NULL) {
    # save default directory info
    defaultDir <- getwd()

    # load and parse sample info
    df.sampleInfo <-
      read.csv(
        file = pathToSampleInfo,
        header = TRUE,
        sep = "\t",
        colClasses = c(
          'character', # outputName
          'character', # htseq_hg38
          'character' # htseq_plasmid
        )
      )
    #return(df.sampleInfo)
    
    
    # load GATC maps
    gatc.hg38 <-
      load.gatc.intervals(intervals.file.path = pathToGatcHg38)
    #return(gatc.hg38)
    
    # TODO: this throws 'Warning: NAs introduced by coercion' due to inability
    # to assign the 'tileID' using the normal 'chrA' parse logic
    #
    # however, the 'tileID' valus isn't used for anything important other than
    # visual sorting confirmation for debugging purposes in developing analyzeRaPID
    #
    # ignore and fix later (or just remove 'tileID' entirely)
    gatc.plasmid <-
      load.gatc.intervals(intervals.file.path = pathToGatcPlasmids)
    #return(gatc.plasmid)
    
    # main sample processing loop
    for (i in 1:nrow(df.sampleInfo)) {
      # set working directory to input
      setwd(dir = inputDataDir)
      
      # parse sample input information
      outputName <- df.sampleInfo[i, 'outputName']
      
      htseq_hg38 <- df.sampleInfo[i, 'htseq_hg38']
      htseq_plasmid <- df.sampleInfo[i, 'htseq_plasmid']
      
      # check files
      print(paste0("hg38 data: ", htseq_hg38))
      print(paste0("plasmid data: ", htseq_plasmid))
      
      # load htseq data
      df.htseq.hg38 <- read.table(
        gzfile(htseq_hg38),
        header = FALSE,
        sep = "\t",
        col.names = c("attribute", "count"),
        stringsAsFactors = FALSE
      )
      
      df.htseq.plasmid <- read.table(
        gzfile(htseq_plasmid),
        header = FALSE,
        sep = "\t",
        col.names = c("attribute", "count"),
        stringsAsFactors = FALSE
      )
      
      
      # calculate basic FKPM for a single samples combining information from
      # reads aligned to hg38 or plasmids
      hg38Processed <- calc.DamID.FKPM.Hg38.withPlasmids(
        data.HTseq.hg38 = df.htseq.hg38,
        gatc.map.hg38 = gatc.hg38,
        data.HTseq.plasmid = df.htseq.plasmid,
        gatc.map.plasmid = gatc.plasmid,
        dsbMask = dskMask
      )
      

      # # make GRanges from data for plasmids
      df.output <-
        select(hg38Processed, seqname, start, end, FKPM)
      
      colnames(df.output) <-
        c("seqname", "start", "end", "score")
      
      gr.rapidCorrected <- makeGRangesFromDataFrame(
        df = df.output,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "seqname",
        start.field = "start",
        end.field = "end",
        seqinfo = Seqinfo(genome = genomeName)
      )
      #return(gr.rapidPlasmid)
      
      
      # set to output directory
      setwd(dir = defaultDir)
      setwd(dir = outputDir)
      
      # write outputs as bigWigs
      rtracklayer::export(
        object = gr.rapidCorrected,
        con = paste0(outputName, ".bw"),
        format = "bigWig"
      )
      
      
      # reset working directory
      setwd(dir = defaultDir)
      
    }
  }

#' 
#' 
#' FOR PLASMID-ALIGNED HT-SEQ PROCESSING
## ----------------------------------------------------------------------------------
#short function to merge GFF & HTSeq data & calculate FKPMs

# REF for FKPM calculation:
# https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

# aggregate data from DamID for reads aligning to hg38 and reads aligning to plasmids

calc.DamID.FKPM.plasmid <-
  function(data.HTseq.hg38,
           gatc.map.hg38,
           data.HTseq.plasmid,
           gatc.map.plasmid,
           dsbMask) {
    # set chromosome order for plasmids
    chr.order <- c(
      "pAC1_HBB_090",
      "pAC2_HBB_080",
      "pAC3_HBB_060",
      "pAC4_HBB_040",
      "pAC5_HBB_020",
      "pAC6_HBB_099",
      "pAC7_HBB_095",
      "pCY29_HBB_100",
      "pCY31_HBB_000",
      "pCY48_HBB_000",
      "pCY49_HBB_100"
    )
    
    # merging of HTSeq with GFF file
    # will exclude the unmapped HTSeq reads (e.g. "__ambiguous", etc.)
    
    data.merged.hg38 <- merge(
      x = data.HTseq.hg38,
      y = gatc.map.hg38,
      by.x = "attribute",
      by.y = "attribute"
    )
    
    data.merged.plasmid <- merge(
      x = data.HTseq.plasmid,
      y = gatc.map.plasmid,
      by.x = "attribute",
      by.y = "attribute"
    )
    
    # ADDED 9 October 2024
    #
    # perform before calculating FKPM to ensure proper denominator rebalancing
    # only perform is a mask is defined
    
    if (!is.null(dsbMask)) {
      # parse 'dsbMask' to (1) seqname; (2) start; (3) end
      dsbParsed = dsbMask %>% strsplit(split = "[:-]") %>% unlist()
      
      
      # select tiles to mask & remove if ANY part of GATC tile overlaps with mask
      maskIndex <- (data.merged.hg38$seqname == dsbParsed[1]) &
        (data.merged.hg38$end >= as.integer(dsbParsed[2])) &
        (data.merged.hg38$start <= as.integer(dsbParsed[3]))
      
      
      # remove masked GATC tiles from dataset and save
      data.merged.hg38 <- data.merged.hg38[!maskIndex, ]
    }
    
    
    # extract totals values (defined at reads to hg38 + reasd to plasmid)
    totalReads <-
      sum(data.merged.hg38$count) + sum(data.merged.plasmid$count)
    #totalReads.plasmid <- data.merged.plasmid$count %>% sum()
    #return(totalReads)
    
    # calculate modified FKPM value as follows:
    # (fragments per 1kbps mapped to plasmid GATC tiles) per (1E6 reads mapped to hg38)
    
    data.merged.plasmid$FKPM <-
      (data.merged.plasmid$count) / ((totalReads / 1E6) * (data.merged.plasmid$length.bp / 1000))
    
    
    
    # formatting output
    data.processed <-
      data.merged.plasmid %>% select("seqname",
                                     "start",
                                     "end",
                                     "tileID",
                                     "length.bp",
                                     "count",
                                     "FKPM")
    
    # seems to output as lexographic for 'tileID' column without this reordering
    data.processed <-
      data.processed[order(match(data.processed$seqname, chr.order),
                           data.processed$tileID), ]
    
    
    return(data.processed)
  }

calculateRapidPlasmidFkpmForPlasmids <-
  function(pathToSampleInfo,
           pathToPlasmidSizes,
           pathToGatcHg38,
           pathToGatcPlasmids,
           inputDataDir,
           outputDir,
           dskMask = NULL) {
    # save default directory info
    defaultDir <- getwd()
    
    # load and parse sample info
    df.sampleInfo <-
      read.csv(
        file = pathToSampleInfo,
        header = TRUE,
        sep = "\t",
        colClasses = c(
          'character', # outputName
          'character', # htseq_hg38
          'character' # htseq_plasmid
        )
      )
    #return(df.sampleInfo)
    
    # build Seqinfo() object for plasmids
    df.plasmidInfo <-
      read.csv(file = pathToPlasmidSizes,
               header = FALSE,
               sep = "\t")
    
    seqinfo.plasmids <-
      Seqinfo(seqnames = df.plasmidInfo[, 1], seqlengths = df.plasmidInfo[, 2])
    #return(seqinfo.plasmids)
    
    
    # load GATC maps
    gatc.hg38 <-
      load.gatc.intervals(intervals.file.path = pathToGatcHg38)
    #return(gatc.hg38)
    
    # TODO: this throws 'Warning: NAs introduced by coercion' due to inability
    # to assign the 'tileID' using the normal 'chrA' parse logic
    #
    # however, the 'tileID' valus isn't used for anything important other than
    # visual sorting confirmation for debugging purposes in developing analyzeRaPID
    #
    # ignore and fix later (or just remove 'tileID' entirely)
    gatc.plasmid <-
      load.gatc.intervals(intervals.file.path = pathToGatcPlasmids)
    #return(gatc.plasmid)
    
    # main sample processing loop
    for (i in 1:nrow(df.sampleInfo)) {
      # set working directory to input
      setwd(dir = inputDataDir)
      
      # parse sample input information
      outputName <- df.sampleInfo[i, 'outputName']
      
      htseq_hg38 <- df.sampleInfo[i, 'htseq_hg38']
      htseq_plasmid <- df.sampleInfo[i, 'htseq_plasmid']
      

      # load htseq data
      df.htseq.hg38 <- read.table(
        gzfile(htseq_hg38),
        header = FALSE,
        sep = "\t",
        col.names = c("attribute", "count"),
        stringsAsFactors = FALSE
      )
      
      df.htseq.plasmid <- read.table(
        gzfile(htseq_plasmid),
        header = FALSE,
        sep = "\t",
        col.names = c("attribute", "count"),
        stringsAsFactors = FALSE
      )
      
      
      # calculate basic FKPM for a single samples combining information from
      # reads aligned to hg38 or plasmids
      plasmidProcessed <- calc.DamID.FKPM.plasmid(
        data.HTseq.hg38 = df.htseq.hg38,
        gatc.map.hg38 = gatc.hg38,
        data.HTseq.plasmid = df.htseq.plasmid,
        gatc.map.plasmid = gatc.plasmid,
        dsbMask = dskMask
      )
      

      # # make GRanges from data for plasmids
      df.output <-
        select(plasmidProcessed, seqname, start, end, FKPM)
      
      
      colnames(df.output) <-
        c("seqname", "start", "end", "score")
      
      gr.rapidPlasmid <- makeGRangesFromDataFrame(
        df = df.output,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "seqname",
        start.field = "start",
        end.field = "end",
        seqinfo = seqinfo.plasmids
      )
      #return(gr.rapidPlasmid)
      
      
      # set to output directory
      setwd(dir = defaultDir)
      setwd(dir = outputDir)
      
      # write outputs as bigWigs
      rtracklayer::export(
        object = gr.rapidPlasmid,
        con = paste0("plasmid_", outputName, ".bw"),
        format = "bigWig"
      )
      
      
      # reset working directory
      setwd(dir = defaultDir)
      
    }
  }

#' 
#' 
#' # STEP 2 - Spike-in correction based on Lambda-gDNA alignment read counts
## ----------------------------------------------------------------------------------
# takes output of calc.DamID.FKPM.plasmid() for 'damid.FKPM.foo'
doSpikeInCorrection <-
  function (bw.rapidFkpm,
            total.counts.signal,
            total.counts.background,
            lambda.counts.signal,
            lambda.counts.background) {
    # calculate lambda factors
    # get lambda reads as % of all reads
    lambda.ratio.signal <-
      (lambda.counts.signal / total.counts.signal)
    
    lambda.ratio.background <-
      (lambda.counts.background / total.counts.background)
    
    
    # ratio of %lambda for (background / signal)
    
    # a HIGHER %lambda corresponds to LESS methylation of sample (gDNA or plasmid)
    # thus, scaling should be such that a LOWER %lambda causes a HIGHER score
    lambda.scaling <-
      (lambda.ratio.background / lambda.ratio.signal)
    
    
    # apply lambda spike-in correction relative to common background
    output <- bw.rapidFkpm
    
    output$score <-
      (bw.rapidFkpm$score) * (lambda.scaling)
    
    return(output)
  }

makeSpikeInCorrectedBigwigs <-
  function(pathToSampleInfo,
           inputDataDir,
           outputDir) {
    # save default directory info
    defaultDir = getwd()
    
    # load and parse sample info
    df.sampleInfo <-
      read.csv(
        file = pathToSampleInfo,
        header = TRUE,
        sep = "\t",
        colClasses = c(
          'character',          # outputName
          'character',          # inputHg38_signal
          'character',          # inputHg38_bkgrd; not used for this step, just for record keeping
          'numeric',          # total_counts_signal
          'numeric',          # total_counts_bkgrd
          'numeric',          # lambda_counts_signal
          'numeric'           # lambda_counts_bkgrd
        )
      )
    #return(df.sampleInfo)
    
    # main sample processing loop
    for (i in 1:nrow(df.sampleInfo)) {
      # set working directory to input
      setwd(dir = inputDataDir)
      
      # parse sample input information
      outputName <- df.sampleInfo[i, 'outputName']
      
      bigwig_signal <- df.sampleInfo[i, 'inputHg38_signal']
      bigwig_background <- df.sampleInfo[i, 'inputHg38_bkgrd']
      
      # name check
      print(paste0("signal file: ", bigwig_signal))
      print(paste0("background file :", bigwig_background))
      
      
      total_counts_signal <-
        df.sampleInfo[i, 'total_counts_signal']
      total_counts_bkgrd <-
        df.sampleInfo[i, 'total_counts_bkgrd']
      
      lambda_counts_signal <-
        df.sampleInfo[i, 'lambda_counts_signal']
      lambda_counts_bkgrd <-
        df.sampleInfo[i, 'lambda_counts_bkgrd']
      
      
      # load FKPM-corrected bigwig data as GRanges
      # GRanges 'score' column is FKPM value
      #
      
      bw.signal <- rtracklayer::import(con = bigwig_signal, format = "bigWig", as = "GRanges")
      
      
      # perform spike-in correction for signal bigWig based
      # on Lambda vs. total read counts
      bw.signal <-
        doSpikeInCorrection(
          bw.rapidFkpm = bw.signal,
          total.counts.signal = total_counts_signal,
          total.counts.background = total_counts_bkgrd,
          lambda.counts.signal = lambda_counts_signal,
          lambda.counts.background = lambda_counts_bkgrd
        )
      
      # # simply load and created a copy for sample tracking purposes
      # pathToBkgrdBigwig <-
      #   paste0(inputDataDir, "/", bigwig_background)
      #
      # bw.bkgrd <-
      #   rtracklayer::import(con = pathToBkgrdBigwig, format = "bigWig", as = "GRanges")
      
      
      # set to output directory
      setwd(dir = defaultDir)
      setwd(dir = outputDir)
      
      # write spikeIn-corrected signal outputs as bigWigs
      rtracklayer::export(
        object = bw.signal,
        con = paste0(outputName, ".bw"),
        format = "bigWig"
      )
      
      
      
      
      # reset working directory
      setwd(dir = defaultDir)
      
    }
  }

#' 
#' 
#' # STEP 3 - background subtraction (+DSB minus noDSB)
## ----------------------------------------------------------------------------------
doBackgroundSubtractionHg38 <-
  function(pathToSampleInfo,
           inputDataDir,
           outputDir) {

    # load and parse sample info
    df.sampleInfo <-
      read.csv(
        file = pathToSampleInfo,
        header = TRUE,
        sep = "\t",
        colClasses = c(
          'character', # outputName
          'character', # bigwig_signal
          'character'  # bigwig_background
        )
        
      )

    # main sample processing loop
    for(i in 1:nrow(df.sampleInfo)){

      # parse sample input information
      outputName <- df.sampleInfo[i, 'outputName']
      
      bigwig_signal <- paste0(inputDataDir,"/",df.sampleInfo[i, 'bigwig_signal'])
      bigwig_background <- paste0(inputDataDir,"/",df.sampleInfo[i, 'bigwig_bkgrd'])

      gr.signal <-
        rtracklayer::import(con = bigwig_signal, format = "bigWig", as = "GRanges")
      
      gr.bkgrd <-
        rtracklayer::import(con = bigwig_background, format = "bigWig", as = "GRanges")

      
      # do background subtraction (+DSB minus no DSB)
      gr.output <- gr.signal %>% granges()
      gr.output$score <- gr.signal$score - gr.bkgrd$score

      
      # write spikeIn-corrected signal outputs as bigWigs
      outputFile <- paste0(outputDir, "/", outputName, ".bw")
      rtracklayer::export(object = gr.output,
                          con = outputFile,
                          format = "bigWig")
      
      print(outputFile)
    }
    
  }

doBackgroundSubtractionPlasmid <-
  function(pathToSampleInfo,
           inputDataDir,
           outputDir) {

    # load 'sampleInfo' definitions
    df.sampleInfo <-
      read.csv(
        file = pathToSampleInfo,
        header = TRUE,
        sep = "\t",
        colClasses = c(
          'character', # outputName
          'character', # bigwig_signal
          'character', # bigwig_bkgrd
          'character'  # plasmid_seqname
        )
        
      )

    
    # main plotting function loop
    for (i in 1:nrow(df.sampleInfo)) {
      # get Path to datafile & starting params

      outputName <- df.sampleInfo[i, 'outputName']
      
      bigwig_signal <- paste0(inputDataDir,"/",df.sampleInfo[i, 'bigwig_signal'])
      bigwig_background <- paste0(inputDataDir,"/",df.sampleInfo[i, 'bigwig_bkgrd'])
      
      # plasmidSeqname <- df.sampleDefinitions[i, 'plasmid_seqname']
      # print(outputName)
      # print(bigwig_signal)
      # print(bigwig_background)
      
      
      # load GRanges object
      gr.plasmidSignal <-
        rtracklayer::import(con = bigwig_signal, format = "bigWig", as = "GRanges")
      
      gr.plasmidBkgrd <-
        rtracklayer::import(con = bigwig_background, format = "bigWig", as = "GRanges")
      
      
      # # select seqname corresponding to 'sampleInfo.txt' defined plasmid
      # gr.plasmidSignal <-
      #   gr.plasmidSignal[seqnames(gr.plasmidSignal) == plasmidSeqname]
      # 
      # gr.plasmidBkgrd <-
      #   gr.plasmidBkgrd[seqnames(gr.plasmidBkgrd) == plasmidSeqname]
      
      print("test2")
      
      # do background subtraction (+DSB minus no DSB)
      gr.plasmidCorrected <- gr.plasmidSignal %>% granges()
      gr.plasmidCorrected$score <-
        gr.plasmidSignal$score - gr.plasmidBkgrd$score
      
      # print maximum GATC tile value for each pair after background subtraction
      print(gr.plasmidCorrected$score %>% max())
      
      
      # write background-subtracted signal outputs as bigWigs
      outputFile <- paste0(outputDir, "/", outputName, ".bw")
      rtracklayer::export(object = gr.plasmidCorrected,
                          con = outputFile,
                          format = "bigWig")
      
      print(outputFile)
    }
    
    
  }

#' 
#' 
#' STEP 3b - bin RaPID-seq to target window size (e.g., 25kbps); ONLY for Hg38
## ----------------------------------------------------------------------------------
doBinRapid <- function(sampleInfo, genomeName, inputDir, outputDir) {
  
  # load 'sampleInfo' definitions
    df.sampleDefinitions <-
      read.csv(
        file = sampleInfo,
        header = TRUE,
        sep = "\t",
        colClasses = c(
          'character', # outputName
          'character', # bigwig_signal, not used here
          'character', # bigwig_bkgrd, not used here
          'character'  # binSizes
        )
      )
    
  # bin samples to defined sizes
for (i in 1:nrow(df.sampleDefinitions)) {
  
  # parse sample pair information
  outputNamePrefix <- df.sampleDefinitions[i, 'outputName']
  
  # recycle from previous step expected filename
  inputFile <-
    paste0(inputDir, "/", df.sampleDefinitions[i, 'outputName'], ".bw")
  
  binSizes <-
    df.sampleDefinitions[i, 'binSizes'] %>% strsplit(";") %>% unlist() %>% as.numeric()
  
  
  # load RaPID-seq raw bigWig file to GRanges

  gr.rapid.raw <-
    rtracklayer::import(con = inputFile, format = "bigWig")
      
      # bin raw DamID differential enrichment values to target resolutions via GRanges::binnedAverage
      
      for (j in 1:length(binSizes)) {
        
        resolution <- binSizes[j]
        outputFile <-
          paste0(outputDir,
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
          gr.RaPID = gr.rapid.raw,
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
    
}

#' 
#' 
#' STEP 4 - make plots
## ----------------------------------------------------------------------------------
# for HG38 aligned datafiles
### single dataplot only
makePlotGeneric <-
  function(df.plot,
           outputPath,
           outputName,
           yMax = 1,
           traceColor = "darkred", linetype = 1) {
    
  # set yMax integer boundary
  yMax <- as.integer(yMax + 1)
  
  # basic params
  xRange <- c(-2E6, +2E6)
  yRange <- c(-1, yMax + 0.5)
  linewidth = 1.5
  
  # customize axis labels
  xlabels <-
    c("-3 Mbps",
      "-2 Mbps",
      "-1 Mbps",
      "0 Mbps",
      "1 Mbps",
      "2 Mbps",
      "3 Mbps")
  
  
  # for h/vline; "gray90" for for solid line; "gray60" for dotted line
  plot <- ggplot(data = df.plot,
                 aes(x = midpoint,
                     y = score)) + geom_hline(
                       yintercept = 0,
                       color = "gray60",
                       linetype = "dotted",
                       linewidth = linewidth
                     ) + geom_vline(
                       xintercept = 0,
                       color = "gray60",
                       linetype = "dotted",
                       linewidth = linewidth
                     ) + geom_area(fill = "#C69180") + geom_line(color = traceColor, linewidth = linewidth, linetype = linetype) + labs(x = "distance from DSB (bps)", y = "RaPID score (raw)") + coord_cartesian(xlim = xRange,
                                                                                                                                                              ylim = yRange,
                                                                                                                                                              expand = TRUE) +  scale_y_continuous(breaks = c(0 , yMax)) +  scale_x_continuous(breaks = c(-3E6, -2E6, -1E6, 0 , 1E6, 2E6, 3E6),
                                                                                                                                                                                                                                             labels = xlabels)
  
  # set further plotting element details
  plot <- plot + theme_classic()
  
  plot <- plot + theme(
    axis.ticks.length = unit(10, "pt"),
    axis.ticks = element_line(size = linewidth, color = "black"),
    axis.line = element_line(linewidth = linewidth, lineend = "square"),
    axis.text = element_text(size = 30,
                             color = "black"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title = element_text(size = 20, color = "black"),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20))
  )

  # save PDF aspect ratio of 4:1
  ggsave(
    filename = paste0(outputPath, "/", outputName, ".pdf"),
    plot = plot,
    device = "pdf",
    width = 6000,
    height = 1500,
    units = "px"
  )
  
  return(plot)
}



### with dotted line for reference data

makePlotGeneric2 <-
  function(df.refData,
           df.plot,
           outputPath,
           outputName,
           yMax = 1,
           traceColorSignal = "darkred",
           traceColorBkgrd = "gray50") {
    # set yMax integer boundary
    yMax <- as.integer(yMax + 1)
    
    # basic params
    xRange <- c(-2E6,+2E6)
    yRange <- c(-1, yMax + 0.5)
    linewidth = 1.5
    
    # customize axis labels
    xlabels <-
      c("-3 Mbps",
        "-2 Mbps",
        "-1 Mbps",
        "0 Mbps",
        "1 Mbps",
        "2 Mbps",
        "3 Mbps")
    
    
    # for h/vline; "gray90" for for solid line; "gray60" for dotted line
    plot <- ggplot(data = df.plot,
                   aes(x = midpoint,
                       y = score)) + geom_hline(
                         yintercept = 0,
                         color = "gray60",
                         linetype = "dotted",
                         linewidth = linewidth
                       ) + geom_vline(
                         xintercept = 0,
                         color = "gray60",
                         linetype = "dotted",
                         linewidth = linewidth
                       ) + geom_area(fill = "#C69180")+ geom_line(color = traceColorSignal, linewidth = linewidth) + labs(x = "distance from DSB (bps)", y = "RaPID score (raw)") + coord_cartesian(xlim = xRange,
                                                                                                                                                                       ylim = yRange,
                                                                                                                                                                       expand = TRUE) +  scale_y_continuous(breaks = c(0 , yMax)) +  scale_x_continuous(breaks = c(-3E6,-2E6,-1E6, 0 , 1E6, 2E6, 3E6),
                                                                                                                                                                                                                                                        labels = xlabels)
    
    # add dotted line of df.refData
    plot <-
      plot + geom_line(
        data = df.refData,
        aes(y = score),
        color = traceColorBkgrd,
        linewidth = linewidth,
        linetype = "11",
        size = 1.5
      )
    
    # set further plotting element details
    plot <- plot + theme_classic()
    
    plot <- plot + theme(
      axis.ticks.length = unit(10, "pt"),
      axis.ticks = element_line(size = linewidth, color = "black"),
      axis.line = element_line(linewidth = linewidth, lineend = "square"),
      axis.text = element_text(size = 30,
                               color = "black"),
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.title = element_text(size = 20, color = "black"),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20))
    )
    
    # save PDF aspect ratio of 4:1
    ggsave(
      filename = paste0(outputPath, "/", outputName, ".pdf"),
      plot = plot,
      device = "pdf",
      width = 6000,
      height = 1500,
      units = "px"
    )
    
    return(plot)
  }


#
#
# for PLASMID aligned datafiles
getRoiTotalScore <- function(gr.rapid, dsbSiteString, radius) {
  # parse DSB site
  dsbParsed <- parseDsbSiteString(dsbSite = dsbSiteString)
  
  
  # setup GRanges for region of interest
  gr.dsb <-
    GRanges(
      seqnames = dsbParsed$chr,
      ranges = IRanges(start = dsbParsed$position - radius, end = dsbParsed$position + radius)
    )
  
  # get overlaps between gr.rapid (genome-wide) and gr.dsb
  gr.overlaps <-
    getDsbOverlaps(gr.data = gr.rapid,
                   gr.dsb = gr.dsb,
                   type = "include") #%>% normalizeRapidLocally
  
  # compute total score
  totalScore <- gr.overlaps$score %>% sum()
  
  
  
  
  # return output
  return(totalScore)
}

#' 
#' 
## ----------------------------------------------------------------------------------
# plotting function for plasmids Rapid score WITHOUT any binning
plotPlasmidRaw <-
  function(gr.plasmidSignal,
           gr.plasmidBkgrd = NULL,
           yMax, 
           seqName) {
    
  # basic params
  linewidth = 1.5
  yRange <- c(-(yMax * 0.067), yMax * 1.067)
  
  # select relevant data
  gr.signal <- gr.plasmidSignal[seqnames(gr.plasmidSignal) == seqName]


  # convert GRanges to dataframe
  gr.signal$midpoint <- gr.signal %>% mid()
  df.signal <- gr.signal %>% as.data.frame()
  
  
  # main plotting function for signal data
  plot <- ggplot(data = df.signal,
                 aes(x = midpoint,
                     y = score)) + geom_hline(
                       yintercept = 0,
                       color = "gray50",
                       linetype = "dotted",
                       linewidth = linewidth
                     ) + geom_area(fill = "#C69180") + geom_line(color = "darkred", linewidth = linewidth) 
  
  
  # add trace corresponding to background (noDSB) if provided
  if (!is.null(gr.plasmidBkgrd)) {
    # select relevant data
    gr.bkgrd <- gr.plasmidBkgrd[seqnames(gr.plasmidBkgrd) == seqName]
    
    # convert GRanges to dataframe
    gr.bkgrd$midpoint <- gr.bkgrd %>% mid()
    df.bkgrd <- gr.bkgrd %>% as.data.frame()  
    
    # add background trace
    plot <-
      plot + geom_line(
        data = df.bkgrd,
        aes(y = score),
        color = "gray30",
        linewidth = linewidth,
        linetype = "11",
      )
  }

  # set theme elements
  plot <- plot + theme_classic()
  
  plot <-
    plot + labs(x = "position (bps)", y = "FKPM") + coord_cartesian(ylim = yRange, expand = TRUE) + scale_y_continuous(breaks = c(0 , yMax))
  
  plot <- plot + theme(
    axis.ticks.length = unit(10, "pt"),
    axis.ticks = element_line(linewidth = linewidth, color = "black"),
    axis.line = element_line(linewidth = linewidth, lineend = "square"),
    axis.text = element_text(size = 30,
                             color = "black"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title = element_text(size = 20, color = "black"),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20))
  ) 
  
  return(plot)
}

makeMaskedGatcMap <- function(pathToGatcBed, tileCutoff = 200) {
  # load GATC map (BED file) to GRanges
  gr.plasmidGatcMap <-
    rtracklayer::import(con = pathToGatcBed, format = "BED")
  
  # set masked bool
  gr.plasmidGatcMap$masked <- (width(gr.plasmidGatcMap) < tileCutoff)
  
  # only keep masked region and merge adjacent
  gr.plasmidGatcMapMasked <-
    gr.plasmidGatcMap[gr.plasmidGatcMap$masked == TRUE] %>% GenomicRanges::reduce()
  
  return(gr.plasmidGatcMapMasked)
}

# Calculate total score-width for RaPID seq on GRanges object 'with or without' 
# applying a mask definition (e.g., GATC tiles < target size)
#
# interpretation equivalent to:
# 'total fragments mapped per kbps of a given region per 1E6 total aligned reads'

calcTotalScore <- function(gr.rapid, gr.gatcMask = NULL){
  
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


#' 
#' 
#' 
