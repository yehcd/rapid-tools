#!/bin/bash

# Charles D Yeh
#
### CHANGES 
#	(5 June 2024)	v0	base version
#	(7 June 2024)	v1	changed 'bamCoveage' to '--normalizeUsing None' for plasmid-samples due to the way filtering affects score
#	(23 Jan 2025)	v2	set bowtie2 genome paths as tunable; cleaned up & remove excessive file generation (bamCoverage, raw alignments).


#
# basic pre-processing of raw Illumina NGS data for RaPID-seq with plasmids.
# ...performs basic trimming of TruSeq read-through, amplicon contaminant filtering, and DamID-adapters
# ...aligns trimmed FASTQs to lambda genome to get spike-in read counts
# ...aligns trimmed FASTQs to plasmid indexes to filter plasmid reads
# 
# ...aligns trimmed+filtered FASTQs to hg38 to get genomic RaPID-seq values
# ...use HTseq to convert to GATC-tile counts FKPM
# 
# ...extract alignment count information for use in analyzeRaPID (R) & other analysis
#
# required packages (conda env with: cutadapt, bowtie2, samtools, deepTools, htseq)
# 
# outputs will be in "./" directory
# runwith "source rapidPreprocessingPlasmid.sh <pathToNgsData> <pathToGatcMap> <inputAmpliconAdr> <pathToPlasmidIndex> <sampleDefinitionFile> <pathToGenomeHg38> <pathToGenomeLambda> <maxThreads>"


printf "\n\nrun RaPID pre-processing on raw NGS data with plasmids\n"

# input argument definitions
# Give as full path to files for proper referencing in conda env! 
# (i.e., use realpath <FILE>)

# $1	directory of NGS FASTQ data
inputNgs=$1

# $2	file: GATC map (hg38)
inputGatcMap=$2

# $3	directory of amplicon filtering sequences
inputAmpliconAdr=$3

# $4	directory containing plasmid bowtie2 index files
inputPlasmidBt2Index=$4

# $5	file: sample root names
sampleDefinitionFile=$5

# $6	path to bowtie2 index genome for hg38
pathToGenomeHg38=$6

# $7	path to bowtie2 index genome for lambda phage
pathToGenomeLambda=$7

# $8	max threads to use
maxThreads=$8

# $9	file: GATC map (plasmid)
plasmidGatcMap=$9




# tunable refernce sequences for cutadapt filtering of amplicon contaminations 
# (R1 = Read1; R2 = Read2)
ampliconSrcR1="../resources/filter_sequences/cutadapt_ampliconPairAdrR1_v2.fa"
ampliconSrcR2="../resources/filter_sequences/cutadapt_ampliconPairAdrR2_v2.fa"




# tunable directories
intermediateFq="fastqs_temp"
outputFastqsKeep="fastqs_amplicons"

outputLambda="bowtie2_Lambda"
outputAmplicon="bowtie2_amplicons"

outputBt2Plasmid="bowtie2_DamID_plasmid"
outputBt2PlasmidCoverage=${outputBt2Plasmid}
outputBt2PlasmidHtseq="${outputBt2Plasmid}/htseq"

outputFilteredHg38="bowtie2_DamID_hg38_plasmidFiltered"
outputFilteredHg38Coverage=${outputFilteredHg38}
outputFilteredHg38Htseq="${outputFilteredHg38}/htseq"



# make directories
printf "\n...making output directories\n"
echo "`date`"

mkdir -p ${intermediateFq} ${outputFastqsKeep} ${outputLambda} ${outputAmplicon} \
${outputBt2Plasmid} ${outputBt2PlasmidCoverage} ${outputBt2PlasmidHtseq} \
${outputFilteredHg38} ${outputFilteredHg38Coverage} ${outputFilteredHg38Htseq}


# main process
{ 
	read
	while IFS="	" read -r output_name fastq_read1 fastq_read2 index_plasmid;
	do 
		printf "\n...processing sample output name: ${output_name}"
		printf "\n...FASTQ Read 1: ${fastq_read1}"
		printf "\n...FASTQ Read 2: ${fastq_read2}\n\n"
		
		# complete path to plasmid bowtie2 index for current sample
		fullPathPlasmidIndex=${inputPlasmidBt2Index}/${index_plasmid}
		
		cutadapt \
		--nextseq-trim=20 \
		-m 20 \
		-j $((maxThreads/4)) \
		--pair-filter=any \
		--times 2 \
		--action=trim \
		-a TruSeq_R1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
		-A TruSeq_R2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
		-o ${intermediateFq}/${output_name}_R1.trim1.fastq.gz \
		-p ${intermediateFq}/${output_name}_R2.trim1.fastq.gz \
		${inputNgs}/${fastq_read1} \
		${inputNgs}/${fastq_read2} \
		> ${outputFastqsKeep}/${output_name}.stderr.cutadapt.trim1_TruseqTrim.log \
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
		-o ${outputFastqsKeep}/${output_name}_R1.amplicons.fastq.gz \
		-p ${outputFastqsKeep}/${output_name}_R2.amplicons.fastq.gz \
		--untrimmed-output ${intermediateFq}/${output_name}_R1.trim2.fastq.gz \
		--untrimmed-paired-output ${intermediateFq}/${output_name}_R2.trim2.fastq.gz \
		${intermediateFq}/${output_name}_R1.trim1.fastq.gz \
		${intermediateFq}/${output_name}_R2.trim1.fastq.gz \
		> ${outputFastqsKeep}/${output_name}.stderr.cutadapt.trim2_ampliconFilter.log \
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
		-o ${intermediateFq}/${output_name}_R1.trim3.fastq.gz \
		-p ${intermediateFq}/${output_name}_R2.trim3.fastq.gz \
		${intermediateFq}/${output_name}_R1.trim2.fastq.gz \
		${intermediateFq}/${output_name}_R2.trim2.fastq.gz \
		> ${outputFastqsKeep}/${i}.stderr.cutadapt.trim3_DamIdTrim.log \
		\
		\
		&& bowtie2 \
		--local \
		--very-sensitive-local \
		--threads ${maxThreads} \
		--phred33 \
		--maxins 800 \
		--soft-clipped-unmapped-tlen \
		-x ${fullPathPlasmidIndex} \
		-1 ${intermediateFq}/${output_name}_R1.trim3.fastq.gz \
		-2 ${intermediateFq}/${output_name}_R2.trim3.fastq.gz \
		2> ${outputBt2Plasmid}/${output_name}.stderr.bt2.DamPlasmid.log\
		\
		| samtools view -S -b - \
		| samtools sort -@ $((maxThreads/4)) - -o ${outputBt2Plasmid}/${output_name}.DamPlasmid.bam \
		\
		&& samtools index ${outputBt2Plasmid}/${output_name}.DamPlasmid.bam \
		&& samtools view -b -f 13 ${outputBt2Plasmid}/${output_name}.DamPlasmid.bam \
		| samtools sort -@ $((maxThreads/4)) -n - -o ${intermediateFq}/${output_name}.unal.DamPlasmid.bam \
		\
		&& bamToFastq -i ${intermediateFq}/${output_name}.unal.DamPlasmid.bam \
		-fq ${intermediateFq}/${output_name}.unal.DamPlasmid.R1.fastq \
		-fq2 ${intermediateFq}/${output_name}.unal.DamPlasmid.R2.fastq \
		\
		&& pigz -p$((maxThreads/4)) ${intermediateFq}/${output_name}.unal.DamPlasmid.R1.fastq \
		&& pigz -p$((maxThreads/4)) ${intermediateFq}/${output_name}.unal.DamPlasmid.R2.fastq \
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
		-1 ${intermediateFq}/${output_name}.unal.DamPlasmid.R1.fastq.gz \
		-2 ${intermediateFq}/${output_name}.unal.DamPlasmid.R2.fastq.gz \
		2> ${outputFilteredHg38}/${output_name}.stderr.bt2.DamHG38.log \
		\
		| samtools view -S -b - \
		| samtools sort -@ $((maxThreads/4)) - -o ${outputFilteredHg38}/${output_name}.DamHG38.bam \
		\
		&& samtools index ${outputFilteredHg38}/${output_name}.DamHG38.bam \
		\
		\
		&& bowtie2 \
		--local \
		--very-sensitive-local \
		--threads ${maxThreads} \
		--phred33 \
		--no-unal \
		--maxins 800 \
		--soft-clipped-unmapped-tlen \
		--no-discordant \
		--no-mixed \
		-x ${pathToGenomeLambda} \
		-1 ${intermediateFq}/${output_name}_R1.trim3.fastq.gz \
		-2 ${intermediateFq}/${output_name}_R2.trim3.fastq.gz \
		2> ${outputLambda}/${output_name}.stderr.bowtie.lambda.log \
		\
		| samtools view -S -b - \
		| samtools sort -@ $((maxThreads/4)) - -o ${outputLambda}/${output_name}.Lambda.bam \
		\
		&& samtools index ${outputLambda}/${output_name}.Lambda.bam \
		\
		&& awk -v name="${output_name}" 'BEGIN {OFS="\t"} /aligned concordantly exactly 1 time/ {print name ".HTseq.out", $1}' "${outputLambda}/${output_name}.stderr.bowtie.lambda.log" \
		>> ./lambda.bt2AlignOnce.counts.txt \
		\
	done } < ${sampleDefinitionFile} 


# HTseq bin-counting en masse for BAM output files
# for plasmid filtered BAM files
for i in `awk '{if(NR>1) print $1}' ${sampleDefinitionFile}`; \
do htseq-count \
--format=bam \
--order=pos \
--stranded=no \
--mode=intersection-nonempty \
--type=region \
--idattr=tileID \
--minaqual=10 \
${outputFilteredHg38}/${i}.DamHG38.bam \
${inputGatcMap} \
>${outputFilteredHg38Htseq}/${i}.HTseq.out 2>${outputFilteredHg38Htseq}/${i}.HTseq.log \
&& gzip ${outputFilteredHg38Htseq}/${i}.HTseq.out \
& done

# wait for above HTseq to finish before proceeding
wait $(jobs -p)


# run HTseq using GATC-map for plasmids
for i in `awk '{if(NR>1) print $1}' ${sampleDefinitionFile}`; \
do htseq-count \
--format=bam \
--order=pos \
--stranded=no \
--mode=intersection-nonempty \
--type=region \
--idattr=tileID \
--minaqual=10 \
${outputBt2Plasmid}/${output_name}.DamPlasmid.bam \
${plasmidGatcMap} \
>${outputBt2PlasmidHtseq}/${i}.HTseq.out 2>${outputBt2PlasmidHtseq}/${i}.HTseq.log \
&& gzip ${outputBt2PlasmidHtseq}/${i}.HTseq.out \
& done

# wait for above HTseq to finish before proceeding
wait $(jobs -p)

# cleanup intermediate datafiles
rm -r ${intermediateFq}
