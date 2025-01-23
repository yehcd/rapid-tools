The directory contains the common functions and datafiles used by the other R Notebook examples for processing of RaPID-seq data. To run these examples, first clone this repository to your local computer and perform the basic preparation steps.


# Preparation
To ensure proper functioning, decompress the following files:

	# Human genome GATC map file (for compatability with 'htseq-count')
	gunzip ./resources/gatcMap_hg38.gff.gz

	# plasmid & lambda phage bowtie2 indexes
	tar -xvzf ./resources/genomes.tar.gz