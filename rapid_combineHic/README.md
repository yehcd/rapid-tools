# RaPID-seq comparisons to chromatin conformation (HiC) data
This example includes the code and sample data for rescaling and comparing RaPID-seq to chromatin conformation data, outputting a plot. This directory is already setup to run; in the following order:

1. combineRapidHic.Rmd
2. makePlots.Rmd


The primary outputs of interest are in:

* `./output/`  


## general data processing workflow
![rapid_hic_workflow](https://github.com/yehcd/rapid-tools/blob/initial/misc/figures/rapid_hic.png)



## description of input datafiles
The datafiles in `./example/` result from processing of raw Illumina NGS data with `../ngs_preprocessing/rapidPreprocessing4.sh` followed by `../rapid_basic/rapidNgsProcessing.Rmd`. The datafiles correspond to:

```
X95_gCY201_12A3_Dox100_binned-25kbps.bw			hg38-aligned, spike-in corrected, background-subtracted, sgHBB

RaPID+HicBalanced_DSBmask-300kbps_individual.RData	pre-combined RaPID/HiC data from all samples as an RData dump
```
