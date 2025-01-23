#!/bin/bash

## Charles Yeh
##
## standard pre-processing & alignment pipeline for RaPID-seq analysis
## removed bowtie2 'min-mapq 10' cutoff as this can be filtered after alignment
## 
## note that in htseq-count, the MAPQ10 cutoff is still set; option '-a 10'
##
## conda environment with the following bioconda packages:
## cutadapt, bowtie2, samtools, deeptools, htseq
##
## inputNames file created using using:
## (cd ngs_data && ls *fastq.gz | sed -e 's/_R[12]_001.fastq.gz//g' | uniq) > names.all.txt
##
## manual adjustment of names files to delete unrelated samples, if necessary
##
## 20 Nov 2023 (v1) - initial version
## 10 Feb 2024 (v2) - added extraction of Lambda count information; added more text and wait for htseq finish
## 15 Feb 2024 (v3) - updated amplicon contamination filtering sequences and added as tunable
## 23 Jan 2025 (v4) - set bowtie2 genome paths as tunable; set maxThreads as tunable



echo "\n\nrun RaPID pre-processing on raw NGS data\n"

# input argument definitions
# Give as full path to files for proper referencing in conda env! 
# (i.e., use realpath <FILE>)

# $1	directory of NGS FASTQ data
inputNgs=$1

# $2	file: GATC map (hg38)
inputGatcMap=$2

# $3	directory of amplicon filtering sequences
inputAmpliconAdr=$3

# $4	file: sample root names
inputNames=$4

# $5	path to bowtie2 index genome for hg38
pathToGenomeHg38=$5

# $6	path to bowtie2 index genome for lambda phage
pathToGenomeLambda=$6

# $7	max threads to use
maxThreads=$7




# tunable refernce sequences for cutadapt filtering of amplicon contaminations 
# (R1 = Read1; R2 = Read2)
ampliconSrcR1="../resources/filter_sequences/cutadapt_ampliconPairAdrR1_v2.fa"
ampliconSrcR2="../resources/filter_sequences/cutadapt_ampliconPairAdrR2_v2.fa"




# tunable directories
outputFQ="fastqs"
outputDam="bowtie2_DamID"
outputLambda="bowtie2_Lambda"
outputAmplicon="bowtie2_amplicons"
outputCoverage=${outputDam}
outputHtseq="${outputDam}/htseq"


# make directories
echo "...making output directories\n"

mkdir -p ${outputFQ}
mkdir -p ${outputDam}
mkdir -p ${outputLambda}
mkdir -p ${outputAmplicon}
mkdir -p ${outputHtseq}


# main alignment processing loop
# trimming for TruSeq, amplicons, DamID followed by alignments
# to hg38 and lambda genome

echo "...start main alignment processing loop"

for i in `cat ${inputNames}`; do 
echo "\n...processing sample name: ${i}\n"
date

cutadapt \
--nextseq-trim=20 \
-m 20 \
-j $((maxThreads/4)) \
--pair-filter=any \
--times 2 \
--action=trim \
-a TruSeq_R1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A TruSeq_R2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o ${outputFQ}/${i}_R1.trimmed.fastq.gz \
-p ${outputFQ}/${i}_R2.trimmed.fastq.gz \
${inputNgs}/${i}_R1*.fastq.gz \
${inputNgs}/${i}_R2*.fastq.gz \
> ${outputFQ}/${i}.stderr.cutadapt.TruseqTrim.log \
\
\
&& cutadapt \
--nextseq-trim=20 \
-m 20 \
-j $((maxThreads/4)) \
-e 0.2 \
--pair-filter=both \
--action=lowercase \
--pair-adapters \
-g "^file:${inputAmpliconAdr}/${ampliconSrcR1}" \
-G "^file:${inputAmpliconAdr}/${ampliconSrcR2}" \
-o ${outputFQ}/${i}_R1.amplicons.fastq.gz \
-p ${outputFQ}/${i}_R2.amplicons.fastq.gz \
--untrimmed-output ${outputFQ}/${i}_R1.ampFiltered.fastq.gz \
--untrimmed-paired-output ${outputFQ}/${i}_R2.ampFiltered.fastq.gz \
${outputFQ}/${i}_R1.trimmed.fastq.gz \
${outputFQ}/${i}_R2.trimmed.fastq.gz \
> ${outputFQ}/${i}.stderr.cutadapt.amplicon.log \
\
\
&& cutadapt \
--nextseq-trim=20 \
-m 20 \
-j $((maxThreads/4)) \
--times 2 \
--pair-filter=any \
--action=trim \
-a R1_DamID_adr_rc=GATCCTCGGCCGCGACC \
-A R2_DamID_adr_rc=GATCCTCGGCCGCGACC \
-g R1_DamID_adr=^GGTCGCGGCCGAGGATC \
-G R2_DamID_adr=^GGTCGCGGCCGAGGATC \
-o ${outputFQ}/${i}_R1.DamID.fastq.gz \
-p ${outputFQ}/${i}_R2.DamID.fastq.gz \
${outputFQ}/${i}_R1.ampFiltered.fastq.gz \
${outputFQ}/${i}_R2.ampFiltered.fastq.gz \
> ${outputFQ}/${i}.stderr.cutadapt.DamID.log \
\
\
&& bowtie2 \
--local \
--very-sensitive-local \
--threads ${maxThreads} \
--phred33 \
--maxins 800 \
--soft-clipped-unmapped-tlen \
--no-discordant \
--no-mixed \
-x ${pathToGenomeHg38} \
-1 ${outputFQ}/${i}_R1.DamID.fastq.gz \
-2 ${outputFQ}/${i}_R2.DamID.fastq.gz \
2> ${outputDam}/${i}.stderr.bowtie2.DamID.log \
\
| samtools view -S -b - \
| samtools sort -@ $((maxThreads/4)) - -o ${outputDam}/${i}.DamID.bam \
\
&& samtools index ${outputDam}/${i}.DamID.bam \
\
\
&& bowtie2 \
--local \
--very-sensitive-local \
--threads ${maxThreads} \
--phred33 \
--no-unal \
--maxins 800 \
--no-discordant \
--no-mixed \
--soft-clipped-unmapped-tlen \
-x ${pathToGenomeLambda} \
-1 ${outputFQ}/${i}_R1.DamID.fastq.gz \
-2 ${outputFQ}/${i}_R2.DamID.fastq.gz \
2> ${outputLambda}/${i}.stderr.bowtie.lambda.log \
\
| samtools view -S -b - \
| samtools sort -@ $((maxThreads/4)) - -o ${outputLambda}/${i}.Lambda.bam \
\
&& samtools index ${outputLambda}/${i}.Lambda.bam \
\
\
&& bamCoverage \
-b ${outputDam}/${i}.DamID.bam \
-p $((maxThreads/4)) \
--normalizeUsing CPM \
--extendReads \
--samFlagInclude 64 \
-o ${outputCoverage}/${i}.DamID.cov.bw \
> ${outputCoverage}/${i}.DamID.cov.stdout+err \
2>&1 \
\
\
&& grep "aligned concordantly exactly 1 time" ${outputLambda}/${i}.stderr.bowtie.lambda.log \
| awk -v name="$i" 'BEGIN {OFS="\t"} {print name ".HTseq.out", $1}' \
>> ./lambda.bt2AlignOnce.counts.txt \
\
\
&& rm ${outputFQ}/${i}_R1.trimmed.fastq.gz \
&& rm ${outputFQ}/${i}_R2.trimmed.fastq.gz \
&& rm ${outputFQ}/${i}_R1.ampFiltered.fastq.gz \
&& rm ${outputFQ}/${i}_R2.ampFiltered.fastq.gz \
&& rm ${outputFQ}/${i}_R1.DamID.fastq.gz \
&& rm ${outputFQ}/${i}_R2.DamID.fastq.gz; \
done

# loop to run htseq jobs in parallel for all samples 
echo "\n...running all HTseq jobs in parallel\n"
date

for i in `cat ${inputNames}`; do
htseq-count \
--format=bam \
--order=pos \
--stranded=no \
--mode=intersection-nonempty \
--type=region \
--idattr=tileID \
--minaqual=10 \
${outputDam}/${i}.DamID.bam \
${inputGatcMap} \
> ${outputHtseq}/${i}.HTseq.out 2> ${outputHtseq}/${i}.HTseq.log \
&& gzip ${outputHtseq}/${i}.HTseq.out & \
done

wait $(jobs -p)
echo "...finished pre-processing!\n"
date
