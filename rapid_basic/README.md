# RaPID-seq comparisons to chromatin conformation (HiC) data
This example includes the code and sample data for converting data from initial raw Illumina NGS processing and performing spike-in correction and background subtraction. Simply run:

```
rapidNgsProcessing.Rmd
```

The primary outputs of interest are in:

* `./output/`  


## general data processing workflow
![rapid_basic_workflow](https://github.com/yehcd/rapid-tools/blob/initial/misc/figures/rapid_basic.PNG)



## description of input datafiles
The datafiles in `./example/` result from processing of raw Illumina NGS data with `../ngs_preprocessing/rapidPreprocessing4.sh`. The datafiles correspond to:

```
319841_46-X104_8_S3.HTseq.out.gz	hg38-aligned, GATC-tile counted, sgHBB
DamID_X104_07_f.HTseq.out.gz		hg38-aligned, GATC-tile counted, sgNTC
```

The spike-in DNA read counts are already embedded in `./example/rapidSampleInfo.txt`.
