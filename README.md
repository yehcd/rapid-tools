# RAD51 proximity identification sequencing (RaPID-seq)
This repository contains the primary code necessary to process next-generation sequencing data from RAD51 proximity identification sequencing (RaPID-seq). 

RaPID-seq is a new method that enables in vivo/in situ monitoring of the homology search process during DNA double-stranded break (DSB) homology-directed repair (HDR). RaPID-seq extends upon the principles of DNA adenine methyltransferase identification (DamID) to provide new information about which HDR donor templates are contacted during repair (candidacy), even if they are not ultimately used as the recombination donor (choice).

These tools generate bigWig datafiles, which can be directly viewed in the [Integrative Genomics Viewer (IGV)](https://igv.org/doc/desktop/).

<br />

![](https://github.com/yehcd/rapid-tools/blob/main/misc/figures/rapid_workflow.png)

<br />

# preparation
To ensure proper functioning, decompress the following files:

```
# Human genome GATC map file (for compatability with 'htseq-count')
gunzip ./resources/gatcMap_hg38.gff.gz

# plasmid & lambda phage bowtie2 indexes
tar -xvzf ./resources/genomes.tar.gz
```

A bowtie2 index for the human genome is also needed (e.g., [UCSC Goldenpath](https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/latest/)).

# usage
The directories in this project are already setup and ordered such that they should directly executable after the above preparation work is completed. In some environments, it may be necessary to use full filepaths to ensure proper functioning. 

