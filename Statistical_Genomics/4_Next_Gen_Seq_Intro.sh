#! /usr/bin/env bash

#use the head command to examine the beginnin of the FASTQ files
##2cells_1.fastq
head 2cells_1.fastq
##2cells_2.fastq
head 2cells_2.fastq
##6h_1.fastq
head 6h_1.fastq
##6h_2.fastq
head 6h_2.fastq

#How many reads are there in each FASQ file?
##2cells_1.fastq
echo $(cat 2cells_1.fastq|wc -l)/4|bc
##2cells_2.fastq
echo $(cat 2cells_2.fastq|wc -l)/4|bc
##6h_1.fastq
echo $(cat 6h_1.fastq|wc -l)/4|bc
##6h_2.fastq
echo $(cat 6h_2.fastq|wc -l)/4|bc

#Run FaastQC for each file
##2cells_1.fastq
~/StatGen/week4/tools/FastQC/fastqc -o ~/StatGen/week4/data/FastQC_output -t 2 --extract --nogroup ~/StatGen/week4/data/2cells_1.fastq
##2cells_2.fastq
~/StatGen/week4/tools/FastQC/fastqc -o ~/StatGen/week4/data/FastQC_output -t 2 --extract --nogroup ~/StatGen/week4/data/2cells_2.fastq
##6h_1.fastq
~/StatGen/week4/tools/FastQC/fastqc -o ~/StatGen/week4/data/FastQC_output -t 2 --extract --nogroup ~/StatGen/week4/data/6h_1.fastq
##6h_2.fastq
~/StatGen/week4/tools/FastQC/fastqc -o ~/StatGen/week4/data/FastQC_output -t 2 --extract --nogroup ~/StatGen/week4/data/6h_2.fastq

#trim the reads in the FASTQ files
##2cells
java -jar ~/StatGen/week4/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 2 -phred33 ~/StatGen/week4/data/2cells_1.fastq ~/StatGen/week4/data/2cells_2.fastq ~/StatGen/week4/data/trimmed/2cells_1.trim.fastq ~/StatGen/week4/data/trimmed/2cells_1.trim.unpaired.fastq ~/StatGen/week4/data/trimmed/2cells_2.trim.fastq ~/StatGen/week4/data/trimmed/2cells_2.trim.unpaired.fastq LEADING:20 TRAILING:20 AVGQUAL:20 MINLEN:25

###Input Read Pairs: 786742 Both Surviving: 770010 (97.87%) Forward Only Surviving: 14831 (1.89%) Reverse Only Surviving: 1596 (0.20%) Dropped: 305 (0.04%)

##6h
java -jar ~/StatGen/week4/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 2 -phred33 ~/StatGen/week4/data/6h_1.fastq ~/StatGen/week4/data/6h_2.fastq ~/StatGen/week4/data/trimmed/6h_1.trim.fastq ~/StatGen/week4/data/trimmed/6h_1.trim.unpaired.fastq ~/StatGen/week4/data/trimmed/6h_2.trim.fastq ~/StatGen/week4/data/trimmed/6h_2.trim.unpaired.fastq LEADING:20 TRAILING:20 AVGQUAL:20 MINLEN:25

###Input Read Pairs: 835648 Both Surviving: 801660 (95.93%) Forward Only Surviving: 29471 (3.53%) Reverse Only Surviving: 3804 (0.46%) Dropped: 713 (0.09%)


#How many reads are there in the trimmed FASTQ files?
##2cells_1.fastq
echo $(cat ~/StatGen/week4/data/trimmed/2cells_1.trim.fastq|wc -l)/4|bc
##2cells_2.fastq
echo $(cat ~/StatGen/week4/data/trimmed/2cells_2.trim.fastq|wc -l)/4|bc
##6h_1.fastq
echo $(cat ~/StatGen/week4/data/trimmed/6h_1.trim.fastq|wc -l)/4|bc
##6h_2.fastq
echo $(cat ~/StatGen/week4/data/trimmed/6h_2.trim.fastq|wc -l)/4|bc

#Count the number of unpaired reads
##2cells_1.fastq
echo $(cat ~/StatGen/week4/data/trimmed/2cells_1.trim.unpaired.fastq|wc -l)/4|bc
##2cells_2.fastq
echo $(cat ~/StatGen/week4/data/trimmed/2cells_2.trim.unpaired.fastq|wc -l)/4|bc
##6h_1.fastq
echo $(cat ~/StatGen/week4/data/trimmed/6h_1.trim.unpaired.fastq|wc -l)/4|bc
##6h_2.fastq
echo $(cat ~/StatGen/week4/data/trimmed/6h_2.trim.unpaired.fastq|wc -l)/4|bc

#Run FaastQC on the trimmed files
##2cells_1.fastq
~/StatGen/week4/tools/FastQC/fastqc -o ~/StatGen/week4/data/FastQC_output -t 2 --extract --nogroup ~/StatGen/week4/data/trimmed/2cells_1.trim.fastq
##2cells_2.fastq
~/StatGen/week4/tools/FastQC/fastqc -o ~/StatGen/week4/data/FastQC_output -t 2 --extract --nogroup ~/StatGen/week4/data/trimmed/2cells_2.trim.fastq
##6h_1.fastq
~/StatGen/week4/tools/FastQC/fastqc -o ~/StatGen/week4/data/FastQC_output -t 2 --extract --nogroup ~/StatGen/week4/data/trimmed/6h_1.trim.fastq
##6h_2.fastq
~/StatGen/week4/tools/FastQC/fastqc -o ~/StatGen/week4/data/FastQC_output -t 2 --extract --nogroup ~/StatGen/week4/data/trimmed/6h_2.trim.fastq


#Index reference genome
~/StatGen/week4/tools/hisat2-2.2.1/hisat2-build ~/StatGen/week4/data/danRer10.chr12.fa danRer10.chr12

#Map the trimmed reads in the 2-cell zebrafish embryo sequencing library
##2cells
~/StatGen/week4/tools/hisat2-2.2.1/hisat2 -q -x ~/StatGen/week4/data/danRer10.chr12 -1 ~/StatGen/week4/data/trimmed/2cells_1.trim.fastq -2 ~/StatGen/week4/data/trimmed/2cells_2.trim.fastq -S ~/StatGen/week4/data/mapped/2cells.sam

##6h
~/StatGen/week4/tools/hisat2-2.2.1/hisat2 -q -x ~/StatGen/week4/data/danRer10.chr12 -1 ~/StatGen/week4/data/trimmed/6h_1.trim.fastq -2 ~/StatGen/week4/data/trimmed/6h_2.trim.fastq -S ~/StatGen/week4/data/mapped/6h.sam


#Convert SAM to BAM files
##2cells
~/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools view -bS ~/StatGen/week4/data/mapped/2cells.sam > ~/StatGen/week4/data/mapped/2cells.bam

##6h
~/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools view -bS ~/StatGen/week4/data/mapped/6h.sam > ~/StatGen/week4/data/mapped/6h.bam

#View the contents of the BAM files
##2cells
~/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools view ~/StatGen/week4/data/mapped/2cells.bam | head

##6h
~/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools view ~/StatGen/week4/data/mapped/6h.bam | head

#Sort and index BAM files
##2cells
~/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools sort -O bam -T tmp_ -o ~/StatGen/week4/data/mapped/2cells_sorted.bam ~/StatGen/week4/data/mapped/2cells.bam
~/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools index ~/StatGen/week4/data/mapped/2cells_sorted.bam

##6h
~/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools sort -O bam -T tmp_ -o ~/StatGen/week4/data/mapped/6h_sorted.bam ~/StatGen/week4/data/mapped/6h.bam
~/StatGen/week4/tools/samtools-bcftools-htslib-1.0_x64-linux/bin/samtools index ~/StatGen/week4/data/mapped/6h_sorted.bam
