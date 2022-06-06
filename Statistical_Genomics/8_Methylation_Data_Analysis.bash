#!/bin/bash

# Data
# Extract data
# This was done in advance, so skip this step
# for INPUT in data/ENCFF857QML_*.fastq.gz; do
# 	echo "extracting file $INPUT"
# 	gzip -dk $INPUT
# done


# 2 Quality control with FastQC
# Create a folder to store the output and run FastQC for each .fastq file
for INPUT in data/ENCFF857QML_*.fastq.gz; do
	~/StatGen/week4/tools/FastQC/fastqc -o results/fastQC -2 --extract --nogroup $INPUT
done

# Interpret plots. Help yourself with FastQC_Manual

# 3 Mapping with Bismark
# 3.1 Preparing the genome
# Download the soft-masked sequence for chromosome 21: Homo_sapiens.GRCh38.dna_sm.chromosome.21.fa.gz and extract it

# The files are consistently named following this pattern:
#    <species>.<assembly>.<sequence type>.<id type>.<id>.fa.gz
# Sequence type:
#   * 'dna_sm' - soft-masked genomic DNA. All repeats and low complexity regions
#     have been replaced with lowercased versions of their nucleic base

echo "Extracting genome file"
gzip -d GRCh38/Homo_sapiens.GRCh38.dna_sm.chromosome.21.fa.gz

# Prepare the genome for mapping
echo "Preparing the genome for mapping"
bismark_genome_preparation --verbose GRCh38/

# What is the content of the .fa files in GRCh38/Bisulfite_genome/CT_conversion and GA_conversion?
# The genome where_
# 	- all the Cs are converted to Ts
# 	- all the Gs are converted to As

# Compute the percentage of “As”, “Cs”, “Gs”, and “Ts” for each of the .fa files.
for i in GRCh38/Bisulfite_Genome/*_conversion/genome_mfa.*_conversion.fa ; do
	echo "Nucleotide content in $i"
	for j in A C G T ; do
		echo -n "$j:"
		grep -o $j $i | wc -l;
	done
done

# Nucleotide content in GRCh38/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa
# A:11820664
# C:1
# G:8226381
# T:20041575
# Nucleotide content in GRCh38/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa
# A:20047046
# C:8185244
# G:1
# T:11856330

# 3.2 Running Bismark
# Process and map the reads using Bismark
echo "Mapping the reads with Bismark"
bismark -N 1 ./GRCh38/ -1 data/ENCFF857QML_R1.fastq.gz -2 data/ENCFF857QML_R2.fastq.gz -o results2/bismark --bam

# Inspect the report file. How do you interpret?

##############################################################################
		# View results/bismark/ENCFF857QML_R1_bismark_bt2_PE_report.txt
##############################################################################

# 4 Extracting methylation calls
# Remove PCR artifacts
echo "Deduplicating"
deduplicate_bismark --bam results2/bismark/ENCFF857QML_R1_bismark_bt2_pe.bam --output_dir results2/bismark/

# Sort and index the dupicated bam diles  using Samtools
	# Sort
echo "Sorting deduplicated file"
samtools sort \
-o results2/bismark/ENCFF857QML_R1_bismark_bt2_pe.sorted.bam \
-O 'bam' \
-T tmp_ results2/bismark/ENCFF857QML_R1_bismark_bt2_pe.deduplicated.bam 

	# Index
echo "Indexing sorted file"
samtools index -b results2/bismark/ENCFF857QML_R1_bismark_bt2_pe.sorted.bam

# View sorted mathylation calls in the bismark file
samtools veiw ENCFF857QML_R1_bismark_bt2_pe.sorted.bam | less

# Extract the methylation call for each "C" in a more readable format
#echo "Extracting methylation call for each C"
bismark_methylation_extractor --report -s \
--counts \
--bedGraph \
--gzip \
-o results2/bismark/Methylation \
results2/bismark/ENCFF857QML_R1_bismark_bt2_pe.sorted.bam

# 4.1 M-bias plot
# The file is generated automatically in the previous step

##############################################################################
# View results/bismark/Methylation/ENCFF857QML_R1_bismark_bt2_pe.sorted.M-bias_R1.png
##############################################################################

# 5 Visualizing methylation data
# 5.1 Generatirng heatmaps and profiles wiht deepTools

# Download the information for the transcription start sites (TSS) in the genome
	# 1) Exclude TSS not corresponding to protein-coding genes and
	# 2) Add a column between TSS and gene stable ID containing the position TSS+1
echo "Processing the TSS file to .bed"
awk '/protein_coding/' results2/bismark/visualization/mart_export.txt | awk -v s=1 '{print $1,$2,$2+s,$3}' OFS='\t' > results2/bismark/visualization/tss.bed

# Generate coverage track
echo "Generating coverage track"
bamCoverage \
--binSize 10 \
--extendReads \
--ignoreDuplicates \
--normalizeUsing RPKM \
--bam results2/bismark/ENCFF857QML_R1_bismark_bt2_pe.sorted.bam \
--outFileFormat bigwig \
--outFileName results2/bismark/visualization/ENCFF857QML_R1_.bigwig \
--numberOfProcessors 2

# Compute matrix for visualization
echo "Computing visualization matrix"
computeMatrix reference-point -S results2/bismark/visualization/ENCFF857QML_R1_.bigwig  \
-R results2/bismark/visualization/tss.bed \
--beforeRegionStartLength=1500 --afterRegionStartLength=1500 \
--referencePoint=center \
-o results2/bismark/visualization/matrix_sample.gz

# Plot heatmap
echo "Plotting heatmap"
plotHeatmap --matrixFile results2/bismark/visualization/matrix_sample.gz \
--outFileName results2/bismark/visualization/heatmap_ENCFF857QML_R1.png

# Plot profile
echo "Plotting profile"
plotProfile --matrixFile results2/bismark/visualization/matrix_sample.gz \
--outFileName results2/bismark/visualization/profile_ENCFF857QML_R1.png
