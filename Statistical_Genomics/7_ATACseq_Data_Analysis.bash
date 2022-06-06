#!/bin/bash

# Preprocessing ATAC-seq data
# 2.1 Indexing BAM-formatted files
# samtools=$(echo "home/up0289/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools")

# Sorting
for INPUT in data/*chr22.bam; do
	OUTPUT=$(echo $INPUT | sed "s/\.bam/\.sorted\.bam/" | sed "s/data\//results\//")
	$samtools sort -o $OUTPUT -O 'bam' -T tmp_ $INPUT
	echo -e "sorted $INPUT file"
done

# Indexing
for INPUT in results/*chr22.sorted.bam; do
	$samtools index -b $INPUT
	echo -e "Indexed $INPUT file"
done

# Count the mapped reads
for INPUT in data/*chr22.bam; do
	echo -e "Number of mapped reads in $INPUT:"
	$samtools view -c -F 260 $INPUT
done
	# Number of mapped reads in results/ENCFF691PRG_chr22.sorted.bam:
	# 2907865
	# Number of mapped reads in results/ENCFF726MIS_chr22.sorted.bam:
	# 3544497
	# Number of mapped reads in results/ENCFF815ABY_chr22.sorted.bam:
	# 3407733

# What does -F 260 mean? How many reads are in the files
# -F FLAG
# Do not output alignments with any bits set in FLAG present in the FLAG field.

# 2.2 Filtering uninformative reads
for INPUT in results/*chr22.sorted.bam; do
	OUTPUT=$(echo $INPUT | sed "s/.sorted\.bam/\.filtered\.bam/")
	$samtools view -q 30 -f 2 -b -h $INPUT > $OUTPUT
	echo -e "Filtered uninformative reads in $INPUT"
done

# How many mapped read pairs are left
for INPUT in results/*chr22.filtered.bam; do
	echo -e "Number of mapped reads left in $INPUT:"
	$samtools view -c -F 260 $INPUT
done
	# Number of mapped reads left in results/ENCFF691PRG_chr22.filtered.bam:
	# 2630010
	# Number of mapped reads left in results/ENCFF726MIS_chr22.filtered.bam:
	# 3228592
	# Number of mapped reads left in results/ENCFF815ABY_chr22.filtered.bam:
	# 3093846

# 2.3 Filtering duplicate reads
for INPUT in results/*chr22.filtered.bam; do
	OUTPUT=$(echo $INPUT | sed "s/.filtered\.bam/\.filtered_duplicates\.bam/")
	METRICS=$(echo $INPUT | sed "s/.filtered\.bam/\_dup_metrics.txt/")
	java -jar picard.jar MarkDuplicates	-REMOVE_DUPLICATES true -I $INPUT -O $OUTPUT -M $METRICS
	echo -e "Filtered duplicate reads in $INPUT file"
done

# How many read pairs are duplicates
for INPUT in results/*chr22.filtered_duplicates\.bam; do
	echo -e "Number of filtered duplicates $INPUT:"
	samtools view -c -F 260 $INPUT
done
	# Number of filtered duplicates results/ENCFF691PRG_chr22.filtered_duplicates.bam:
	# 2357280
	# Number of filtered duplicates results/ENCFF726MIS_chr22.filtered_duplicates.bam:
	# 2902316
	# Number of filtered duplicates results/ENCFF815ABY_chr22.filtered_duplicates.bam:
	# 2725826


# 2.4 Checking the size of the fragments
for INPUT in results/*chr22.filtered_duplicates\.bam; do
	OUTPUT=$(echo $INPUT | sed "s/.filtered_duplicates\.bam/\.insert_size_metrics\.txt/")
	HIST=$(echo $INPUT | sed "s/.filtered_duplicates\.bam/\.insert_size_histogram\.pdf/")
	java -jar picard.jar CollectInsertSizeMetrics -I $INPUT -O $OUTPUT -H $HIST
done

############################################################
				# Look at generated pdfs
############################################################

# 3 Peak calling
# Identify the regions where reads are statistically enriched
# macs2=$(echo "/home/up0289/.local/bin/macs2")
for INPUT in results/*chr22.filtered_duplicates\.bam; do
	NAME=$(echo $INPUT | sed "s/.filtered_duplicates\.bam/./" | sed "s/results\/./")
	$macs2 callpeak -t $INPUT \
	-n $NAME \
	--outdir results/enriched_regions/ \
	-f BAMPE \
	-g 2.9e9 \ 
	--nomodel \
	--call-summits
done
# -g mappable genome size - has to do with repeptitive regions (size of the genome without the huge repetitive sequences)

# How many peaks do you obtain from each replicate? Visualize with R
for INPUT in results/enriched_regions/*chr22_peaks\.xls; do
	echo -e "Peaks in $INPUT"
	grep "fragments after filtering in treatment:" $INPUT
done 
	# Peaks in results/enriched_regions/ENCFF691PRG_chr22_peaks.xls
	# fragments after filtering in treatment: 1178640
	# Peaks in results/enriched_regions/ENCFF726MIS_chr22_peaks.xls
	# fragments after filtering in treatment: 1451158
	# Peaks in results/enriched_regions/ENCFF815ABY_chr22_peaks.xls
	# fragments after filtering in treatment: 1362913
	
############################################################
					# Visualize with R
############################################################

# Peaks that overlap across replicates
	# homer=$PATH:/home/up0289/progs/homer/bin/
	# mergePeaks=/home/up0289/progs/homer/bin/mergePeaks
for DIST in 1 25 50 100
do
	OUTPUT=$(echo "results/overlap_peaks/MergedPeakFile_${DIST}.txt")
	~/progs/homer/bin/mergePeaks -d $DIST results/enriched_regions/ENCFF691PRG_chr22_peaks.narrowPeak results/enriched_regions/ENCFF726MIS_chr22_peaks.narrowPeak results/enriched_regions/ENCFF815ABY_chr22_peaks.narrowPeak > $OUTPUT
done

# Count the number of resulting peaks
for INPUT in results/overlap_peaks/MergedPeakFile_*.txt; do
	wc -l $INPUT
done
	# 7002 results/overlap_peaks/MergedPeakFile_1.txt
	# 4465 results/overlap_peaks/MergedPeakFile_100.txt
	# 5972 results/overlap_peaks/MergedPeakFile_25.txt
	# 5277 results/overlap_peaks/MergedPeakFile_50.txt

# Generate two bed files from the merged file for d=50
tail -n+2 results/overlap_peaks/MergedPeakFile_50.txt | awk -v OFS='\t' ' NF==9 {print $2,$3,$4,$1} ' > results/overlap_peaks/MergedPeakFile_50_in_one.bed
tail -n+2 results/overlap_peaks/MergedPeakFile_50.txt | awk -v OFS='\t' ' NF>9 {print $2,$3,$4,$1} ' > results/overlap_peaks/MergedPeakFile_50_in_multiple.bed


# How many peaks do the each of the resulting consensus peak sets consist of
for INPUT in results/overlap_peaks/MergedPeakFile_50_in_*.bed;do
	wc -l $INPUT
done
	# 1297 results/overlap_peaks/MergedPeakFile_50_in_multiple.bed
	# 3979 results/overlap_peaks/MergedPeakFile_50_in_one.bed


# Visualization with deepTools
for INPUT in results/overlap_peaks/MergedPeakFile_50_in_*.bed; do
	$samtools index -b $INPUT
	echo -e "Indexed $INPUT file"
done

for INPUT in results/*chr22.filtered_duplicates.bam; do
	NAME=$(echo $INPUT | sed "s/.filtered_duplicates\.bam/.bigwig/")
	bamCoverage \
	--binSize 10 \
	--extendReads \
	--ignoreDuplicates \
	--normalizeUsing RPKM \
	--bam $INPUT \
	--outFileFormat bigwig \
	--outFileName $NAME \
	--numberOfProcessors 2
done


# Matrix for visualization
for INPUT in results/overlap_peaks/MergedPeakFile_50_in_*.bed; do
	OUTPUT=$(echo $INPUT | sed "s/MergedPeakFile_/vis_matrix_d/" | sed "s/overlap_peaks\//vis_deeptools\//" | sed "s/\.bed/\.gz/")
	echo -e "Computing matrix for $OUTPUT"
	computeMatrix reference-point -S results/ENCFF{691PRG,726MIS,815ABY}_chr22.bigwig \
	-R $INPUT \
	--beforeRegionStartLength=1500 --afterRegionStartLength=1500 \
	--referencePoint=center \
	-o $OUTPUT
done


# Plot heat map
for INPUT in results/vis_deeptools/vis_matrix_d50_in_*.gz; do
	OUTPUT=$(echo $INPUT | sed "s/vis_matrix/heatmap/" | sed "s/\.gz/\.png/")
	plotHeatmap --matrixFile $INPUT --heatmapWidth 10 --outFileName $OUTPUT
done

############################################################
					# Check heatmaps
############################################################

# Motif analysis
cd ~/progs/homer/bin
for INPUT in ~/StatGen/lab7/results/overlap_peaks/MergedPeakFile_50_in_*.txt; do
	./findMotifsGenome.pl $INPUT hg38 ~/StatGen/lab7/results/motif_analysis -size 200 -mask
done


# Explore generated html file