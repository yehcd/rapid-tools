# RaPID-seq with plasmids
This example includes the code and sample data for running the RaPID-seq with plasmids analysis and outputting both bigWig and plot files. This directory is already setup to run; in the following order:

1. dataPreparation.Rmd
2. plot_hg38.Rmd


The primary outputs of interest are in:

* `./output/output3_bkgrdSubtracted`  
* `./output/output4_plots`



## description of input datafiles
The datafiles in `./example/` result from processing of raw Illumina NGS data with `../ngs_preprocessing/rapidPreprocessingPlasmid2.sh`. The datafiles correspond to:

```
C3R5X126_A1_gCY321_noPlasmid.HTseq.out.gz		hg38-aligned, no plasmid, sgNTC
C3R5X126_B1_gCY201_noPlasmid.HTseq.out.gz		hg38-aligned, no plasmid, sgHBB

C3R5X126_A2_gCY321_pCY29.HTseq.out.gz			hg38-aligned, HBB 99.9% plasmid, sgNTC
C3R5X126_B2_gCY201_pCY29.HTseq.out.gz			hg38-aligned, HBB 99.9% plasmid, sgHBB

plasmid_C3R5X126_A1_gCY321_noPlasmid.HTseq.out.gz	plasmid-aligned, no plasmid, sgNTC
plasmid_C3R5X126_B1_gCY201_noPlasmid.HTseq.out.gz	plasmid-aligned, no plasmid, sgHBB

plasmid_C3R5X126_A2_gCY321_pCY29.HTseq.out.gz		plasmid-aligned, HBB 99.9% plasmid, sgNTC
plasmid_C3R5X126_B2_gCY201_pCY29.HTseq.out.gz		plasmid-aligned, HBB 99.9% plasmid, sgHBB
```
