# initial NGS raw data processing

These are simple bash scripts for automated pre-processing of raw FASTQ data from Illumina sequencing. The main steps are:

* cutadapt removal Truseq sequences, amplicon contimants, and DamID sequences
* bowtie2 alignment to human genome (hg38)
* bowtie2 alignment to spike-in control (lambda phage)
* htseq-count of reads into GATC tiles


## requirements

conda environment with the following packages and dependencies:
	cutadapt
	bowtie2
	samtools
	deeptools
	htseq

The bowtie2 index for humans (e.g., hg38) needs to downloaded and prepared first (e.g., [UCSC Goldenpath](https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/latest/)).



## example for standard RaPID-seq NGS pre-processing
```
source ./rapidPreprocessing4.sh \
./example/ \
../resources/gatcMap_hg38.gff \
../resources/filter_sequences/ \
./example/sampleInfo.txt \
\
<HUMAN BOWTIE2 INDEX> \
\
../resources/genomes/lambda/ncbi_Lambda \
60
```

![ngs_standard_workflow](https://github.com/yehcd/rapid-tools/blob/initial/misc/figures/preprocessing_basic.PNG)


## example for standard RaPID-seq NGS pre-processing
For RaPID-seq experiments with plasmids, the plasmid reads are first extracted before performing alignments to the human genome.

```
source ./rapidPreprocessingPlasmid2.sh \
./example/ \
../resources/gatcMap_hg38.gff \
../resources/filter_sequences/ \
../resources/genomes/plasmid/all_indexes \
./example/sampleInfo.txt \
\
<HUMAN BOWTIE2 INDEX> \
\
../resources/genomes/lambda/ncbi_Lambda \
60
../resources/gatcMap_plasmids.gff
```

![ngs_plasmid_workflow](https://github.com/yehcd/rapid-tools/blob/initial/misc/figures/preprocessing_plasmid.PNG)
