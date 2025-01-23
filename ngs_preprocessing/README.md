These are simple bash scripts for automated pre-processing of raw FASTQ data from Illumina sequencing. The main steps are:

* cutadapt removal Truseq sequences, amplicon contimants, and DamID sequences
* bowtie2 alignment to human genome (hg38)
* bowtie2 alignment to spike-in control (lambda phage)
* htseq-count of reads into GATC tiles


Requirements:

conda environment with the following packages and dependencies:
	cutadapt
	bowtie2
	samtools
	deeptools
	htseq

bowtie2 index for humans (e.g., hg38)



# example command for stamdard RaPID-seq alignment of FASTQ data
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


# example command for plasmid RaPID-seq alignment of FASTQ data
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


