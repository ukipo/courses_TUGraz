#!/bin/bash

# 2 Quality control
## Create folder to store FastQC output
mkdir -p fastqc

## Run FastQC for each .fastq file
for INPUT in data/HiC_S2_1p_10min_lowU_R*.fastq.gz; do
	fastqc -o fastqc -t 2 $INPUT
done

############################################################
				# Examine .html files
############################################################

# 3 Mapping the Hi-C reads to the reference genome
# 3.1 Indexing the genome
## Download the doft-masked top-level sequence file of the BDGP6 assembly of the Drosophila melanogaster

## Create a folder called BDGP6 and store the uncompressed file in this folder
mkdir -p BDGP6

## Index the genome for mapping
bowtie2-build BDGP6/Drosophila_melanogaster.BDGP6.32.dna.chromosome.fa.gz \
	./BDGP6/Drosophila_melanogaster.BDGP6.32.dna.chromosome \
	--threads 2


# 3.2 Running Bowtie2
## Generate a directory Mappings
mkdir -p Mappings

## Map the forward and reverse reads
bowtie2 -x ./BDGP6/Drosophila_melanogaster.BDGP6.32.dna.chromosome \
		-U data/HiC_S2_1p_10min_lowU_R1.fastq.gz --threads 2 \
		--reorder --very-sensitive-local | samtools view -Shb \
		> Mappings/HiC_S2_1p_10min_lowU_R1.bam

## Do the same for the reverse reads
bowtie2 -x ./BDGP6/Drosophila_melanogaster.BDGP6.32.dna.chromosome \
		-U data/HiC_S2_1p_10min_lowU_R2.fastq.gz --threads 2 \
		--reorder --very-sensitive-local | samtools view -Shb \
		> Mappings/HiC_S2_1p_10min_lowU_R2.bam

## In a loop
for INPUT in data/HiC_S2_1p_10min_lowU_R*.fastq.gz; do
	OUTPUT=$(echo $INPUT | sed "s/\.fastq.\gz/\.bam/" | sed "s/data\//Mappings\//")
	bowtie2 -x ./BDGP6/Drosophila_melanogaster.BDGP6.32.dna.chromosome \
		-U $INPUT --threads 2 \
		--reorder --very-sensitive-local | samtools view -Shb \
		> $OUTPUT

	## What is the purpose of the --reorder parameter?
	## Output SAM are printed in order corresponding to the order of reads in original input file
	
	## Why local instead of end-to-end?
	## local does nott need the entire alignement from end to end
	
	## Difference between --local and --very-sensitive-local? Why the latter?
	## The latter is slower, but has increased likelyhood of reporting a correct alignement
	## --local: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
	## --very-sensitive-local: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75
	## -D limit seed extensions
	## -R max time it re-seeds with repetitive reads
	## -N amount of mismatches permitted per seed
	## -L length of seed subsetting
	## -i interval between seed subsetting to use in multiseed alignement
	
	## Purpose of Samtools
	## convert to BAM
	## -S autodetects input format
	## -h includes header in SAM
	## -b output BAM
	
	## Percent of successfully mapped reads?


# 4 Constructing the Hi-C contact matrix
## DpnII restriction site sequence?
## GATC

## Sequence left at the DpnII restriction site?
## GATC

## Find restriction sites in the reference genome
hicFindRestSite --fasta BDGP6/Drosophila_melanogaster.BDGP6.32.dna.chromosome.fa.gz \
	--searchPattern GATC -o BDGP6/restriction_sites.bed

## Create folder store matrix of interactions
mkdir -p hicMatrix

## Create matrix of interactions
hicBuildMatrix --samFiles Mappings/HiC_S2_1p_10min_lowU_R1.bam \
	Mappings/HiC_S2_1p_10min_lowU_R2.bam \
	--binSize 10000 \
	--restrictionSequence GATC \
	--restrictionCutFile BDGP6/restriction_sites.bed \
	--danglingSequence GATC \
	--outBam hicMatrix/HiC_S2_1p_10min_lowU_valid.bam \
	--outFileName hicMatrix/HiC_S2_1p_10min_lowU_10kbp.h5 \
	--QCfolder hicMatrix/ \
	--threads 2


############################################################
				# Examine .html file in QC folder
############################################################

	## How many read pairs are considered valid?
	## 5443048	

############################################################
				# Verify valid reads in IGV
############################################################

# 4.1 Visualize the Hi-C matrix
## Merge matrix bins
	## how many bins N to go from 10kbp to 100kbp resolution?
	## 10

hicMergeMatrixBins \
	--matrix hicMatrix/HiC_S2_1p_10min_lowU_10kbp.h5 \
	--numBins 10 \
	--outFileName hicMatrix/HiC_S2_1p_10min_lowU_100kbp.h5


## Generate directory to visualize the merged matrix
mkdir -p Plots

## Visualize the merged matrix
hicPlotMatrix \
	--matrix hicMatrix/HiC_S2_1p_10min_lowU_100kbp.h5 \
	--log \
	--dpi 300 \
	--clearMaskedBins \
	--chromosomeOrder 2L 2R 3L 3R 4 X Y \
	--colorMap jet \
	--title "Hi-C contact map 100kbp" \
	--outFileName Plots/HiC_S2_1p_10min_lowU_100kbp.png


############################################################
				# View the contact maps
############################################################

	## What is plotted on the x and y axis?
	## The positions in the genome, separated by chromosome
	
	## What does a high value far off from the diagonal imply?
	## Interaction between the two regions
	

# 5 Correcting the Hi-C matrix
## Generate histogram
hicCorrectMatrix diagnostic_plot \
	--chromosomes 2L 2R 3L 3R 4 X Y \
	--matrix hicMatrix/HiC_S2_1p_10min_lowU_100kbp.h5 \
	--plotName Plots/HiC_S2_1p_10min_lowU_100kbp.h5_diagnostic_plot.png


############################################################
				# Examine histogram and decide tresholds
############################################################

	## INFO:hicexplorer.hicCorrectMatrix:Removing 4 zero value bins
	## INFO:hicexplorer.hicCorrectMatrix:mad threshold -2.4199773582474227


## Remove bins that do not satisfy your chosen treshold
hicCorrectMatrix correct \
	--matrix hicMatrix/HiC_S2_1p_10min_lowU_100kbp.h5 \
	--chromosomes 2L 2R 3L 3R X \
	--perchr \
	--correctionMethod ICE \
	--filterThreshold -1.6 1.8 \
	--outFileName hicMatrix/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4.h5


## Visualize corrected matrix
hicPlotMatrix \
	--matrix hicMatrix/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4.h5 \
	--log \
	--dpi 300 \
	--clearMaskedBins \
	--chromosomeOrder 2L 2R 3L 3R X \
	--colorMap jet \
	--title "Corrected Hi-C contact map 100kbp" \
	--outFileName Plots/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4.png


# 6 Determining A and B compartments
## Perform this process of 1) generating the pearson correaltion matrix and then 2) extracting the first pronciple component with the function hicPCA to identify A and B compartments in our HiC data set.

## Run hicPCA adjusting the parameters:
hicPCA --matrix hicMatrix/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4.h5 \
	--whichEigenvectors "1" \
	--pearsonMatrix hicMatrix/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4_pearson.h5 \
	--format bigwig \
	--ignoreMaskedBins \
	--outputFileName hicMatrix/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4_pc1.bw


## Visualize the Pearson matrix with the PC scores
hicPlotMatrix --matrix hicMatrix/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4_pearson.h5 \
	--perChromosome \
	--chromosomeOrder 2L 2R 3L 3R X \
	--colorMap hot \
	--bigwig hicMatrix/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4_pc1.bw \
	--title "Pearson matrix and PC1" \
	--outFileName Plots/HiC_S2_1p_10min_lowU_100kbp_corrected-no-y4_pca1.png

	## What does PC1 represent?
	## The maps are separated by chromosomes, to get rid of inter-chromosome interactions
	

# 7 Calling topologically associating domains (TADs)
## Call TADs for the higher resolution matrix (10kbp)

############################################################
# First you have to correct the higher resolution matrix!
############################################################

	## Remove bins that do not satisfy your chosen treshold
	hicCorrectMatrix correct \
		--matrix hicMatrix/HiC_S2_1p_10min_lowU_10kbp.h5 \
		--chromosomes 2L 2R 3L 3R X \
		--perchr \
		--correctionMethod ICE \
		--filterThreshold -1.6 1.8 \
		--outFileName hicMatrix/HiC_S2_1p_10min_lowU_10kbp_corrected.h5


	## Visualize corrected matrix
	#hicPlotMatrix \
	#	--matrix hicMatrix/HiC_S2_1p_10min_lowU_10kbp_corrected.h5 \
	#	--log \
	#	--dpi 300 \
	#	--clearMaskedBins \
	#	--chromosomeOrder 2L 2R 3L 3R 4 X Y \
	#	--colorMap jet \
	#	--title "Corrected Hi-C contact map 10kbp" \
	#	--outFileName Plots/HiC_S2_1p_10min_lowU_10kbp_corrected.png
	
	## Killed <- too big, that's why we do 100kbp above!
	

## Call TADs for the higher resolution matrix (10kbp)
hicFindTADs --matrix ./hicMatrix/HiC_S2_1p_10min_lowU_10kbp_corrected.h5 \
	--minDepth 30000 \
	--maxDepth 100000 \
	--step 10000 \
	--correctForMultipleTesting fdr \
	-- thresholdComparisons 0.05 \
	--delta 0.001 \
	--outPrefix ./hicMatrix/HiC_S2_1p_10min_lowU_10kbp_corrected_TADs \
	--numberOfProcessors 2


## Visualize the TADs of a region of your choice
hicPlotTADs --tracks tracks.ini -o <output_file_name> \
	--region <chr_region> --height 6

## make track files
make_tracks_file \
	--trackFiles hicMatrix/HiC_S2_1p_10min_lowU_10kbp_corrected.h5 \
	./findTADs/HiC_S2_1p_10min_lowU_10kbp_corrected_TADs_tad_score.bm \
	./hicMatrix/HiC_S2_1p_10min_lowU_10kbp_corrected_pc1.bw \
	--out ./hicMatrix/tracks_close_to_tutorial.ini


## Plot
hicPlotTADs --tracks hicMatrix/tracks1.ini \
	-o Plots/HiC_S2_1p_10min_lowU_10kbp_corrected_TADs_all-at-once.png \
	--region 3L:7000000-12000000 \
	--title "Hi-C TADs domains: region 3L:7Mbp-12Mbp" 
	
hicPlotTADs --tracks hicMatrix/tracks1.ini \
	-o Plots/HiC_S2_1p_10min_lowU_10kbp_corrected_TADs_2L_145-165.png \
	--region 2L:14500000-16500000 \
	--title "Hi-C TADs domains: region 2L:145Mbp-165Mbp" 

hicPlotTADs --tracks hicMatrix/tracks1.ini \
	-o Plots/HiC_S2_1p_10min_lowU_10kbp_corrected_TADs_2L_0-5.png \
	--region 2L:00000-500000 \
	--title "Hi-C TADs domains: region 2L:0Mbp-5Mbp" 

hicPlotTADs --tracks hicMatrix/tracks1.ini \
	-o Plots/HiC_S2_1p_10min_lowU_10kbp_corrected_TADs_X_11-15.png \
	--region X:1200000-1500000 \
	--title "Hi-C TADs domains: region X:12Mbp-15Mbp" 











