#!/bin/bash

# Goals
    ## Identify larger structural variants in three individuals sequenced in The 1000 Genomes Project. The subjects are two parents and a child.

# 2 Mapping with BWA-MEM
    ## 1. Download the bam files from the teach center
    ## The reads were already mapped

    ## bwa mem <REF_ASSEMBLY.fa> \
    ##   Reads_Forward.fq \
    ##   Reads_Reverse.fq \
    ##   -M -t 2 | samtools view -S -b -

    ## 2. Download the reference assembly
        ## 1) connect to the ensembl server:
        ftp ftp.ensemble.org

        ## 2) navigate to the correct directory with the cd command
        cd pub/release-102/fasta/homo_sapiens/dna/

        ## 3) download the assembly with the get command
        get Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa.gz

    ## 3. Index the reference sequence
        gunzip -k refgenome/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa.gz
        samtools faidx refgenome/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa
    
    ## 4. Use samtools to sort and index the bam files
    ## Sort
        for INPUT in data/NA128*_chr20.bam; do
            OUTPUT=$(echo $INPUT | sed "s/\.bam/\.sorted\.bam/")
            samtools sort -O 'bam' -T tmp_ -o $OUTPUT $INPUT
        done

    ## Index
        for INPUT in data/NA128*_chr20.sorted.bam; do
            samtools index -b $INPUT
        done

# 3 Characterizing the fragment size distribution
    ## Characterize the insert size distribution. We will use it to decide on the treshold for discordant mappings.

    ## 5. Using picard: The following commands extract pairs of forward and reverse reads from a BAM file and compute the mean and standard deviation of the insert size.
        for INPUT in data/NA128*_chr20.sorted.bam; do
            OUTPUT=$(echo $INPUT | sed "s/\.bam/\.insert_size_metrics\.txt/" | sed "s/data\//frag_size_dist\//")
            HISTOGRAM=$(echo $INPUT | sed "s/\.bam/\.insert_size_histogram\.pdf/" | sed "s/data\//frag_size_dist\//")
            picard CollectInsertSizeMetrics \
                I=$INPUT \
                O=$OUTPUT \
                HISTOGRAM_FILE=$HISTOGRAM
        done
        
        ## Picard tool is pretty nice fro this

        ## a) Use R to visualize the insert size distribution. Repeat for all remaining samples.
            
            ## USE ABSOLUTE VALUE

            R
                library(ggplot2)

                myFiles <- Sys.glob("frag_size_dist/NA128*_chr20.sorted.insert_size_metrics.txt")

                for(i in 1:length(myFiles)){
                    # read in table
                    fragsize <- read.table(file = myFiles[i], header = TRUE, skip = 10)
                    # output
                    outFile <- gsub("NA128", "frag_size_dist/NA128", gsub("metrics.txt", "histogram_R.png", basename(myFiles[i])))
                    # sample name (for title)
                    sampleName <- gsub(".sorted.insert_size_metrics.txt", "", basename(myFiles[i]))
                    # write file
                    png(filename = outFile)
                        # plot
                        print(ggplot(data = fragsize, mapping = aes(x=insert_size, All_Reads.fr_count)) +
                            geom_bar(stat = "identity", color="slateblue1", fill="slateblue1") +
                            labs(title=sampleName))
                    dev.off()
                    # mean and sd of insert sizes
                    cat(c(sampleName, "mean insert size:", mean(fragsize$All_Reads.fr_count), "\n", sampleName, "SD of insert size:", sd(fragsize$All_Reads.fr_count)))
                    }
            q()

        ## Do ABSOLUTE VALUES and trim or you get skewed averages and sd!!!!!!!!!!!!!!!!!!!

        ## b) What is the mean and standard deviation of the insert size? Is it the same for all samples?
            ## NA12878_chr20 mean insert size:  13505.0377604167 
            ## NA12878_chr20 SD of insert size: 19732.7137341544

            ## NA12891_chr20 mean insert size: 14729.5965665236 
            ## NA12891_chr20 SD of insert size: 22147.8405175446

            ## NA12892_chr20 mean insert size: 15581.3789173789 
            ## NA12892_chr20 SD of insert size: 23334.9749713752

    ## 6. Using samtools: The ninth column of a SAM file, observed Template LENgth (TLEN), can be used as an approximate of the fragment length. Let us obtain the insert sizes using only the first pair of properly mapped pairs (flag -f66, see Decoding SAM flags).
        for INPUT in data/NA128*_chr20.sorted.bam; do
            OUTPUT=$(echo $INPUT | sed "s/\.bam/\.insert_sizes\.txt/" | sed "s/data\//frag_size_dist\//")
            samtools view -f66 $INPUT | cut -f 9 > $OUTPUT
        done
        
        samtools view -f66 input.bam | cut -f 9 > insert-sizes.txt

        ## a) Why are there negative values?
            ## TLEN: the number of bases covered by the reads from the same fragment. Plus/minus means the current read is the leftmost/rightmost read.

        ## b) Use R to visualize the distribution of insert sizes and get some summary statistics.
            R
                library(ggplot2)

                myFiles <- Sys.glob("frag_size_dist/NA128*_chr20.sorted.insert_sizes.txt")

                for(i in 1:length(myFiles)){
                    # read in table
                    insertSizes <- read.table(file = myFiles[i])
                    # output
                    outFile <- gsub("NA128", "frag_size_dist/NA128", gsub(".txt", "_distribution.png", basename(myFiles[i])))
                    # sample name (for title)
                    sampleName <- gsub(".sorted.insert_sizes.txt", "", basename(myFiles[i]))
                    # write file
                    png(filename = outFile)
                        # plot
                        print(ggplot() +
                                geom_boxplot(aes(y=insertSizes$V1)) +
                                coord_cartesian(ylim = c(-1500, 1500)) +
                                labs(title=sampleName)
                            )
                    dev.off()
                    # summaries
                    cat(c(sampleName, "summary statistics:", summary(insertSizes)))
                }
            q()

            ## NA12878_chr20 summary statistics:
                ## Min.   :-63455993
                ## 1st Qu.:     -319   
                ## Median :      -50   
                ## Mean   :     -172   
                ## 3rd Qu.:      319   
                ## Max.   : 62505943 

            ## NA12891_chr20 summary statistics: 
                ## Min.   :-59266645   
                ## 1st Qu.:     -300   
                ## Median :      -76   
                ## Mean   :     -163   
                ## 3rd Qu.:      299   
                ## Max.   : 60388879  

            ## NA12892_chr20 summary statistics: 
                ## Min.   :-60530296   
                ## 1st Qu.:     -303   
                ## Median :      -78   
                ## Mean   :     -385   
                ## 3rd Qu.:      302   
                ## Max.   : 61194790  

    ## 7. The insert size estimation using this method has a limitation: it merely reflects the distance between the mappings. In particular, the method may provide misleading estimates for RNA-seq data (Why?).
        ## RNA-seq sequences processed mRNA, which has only exons and no introns or other nc regions.

    ## 8. Exclude very large, unlikely inserts and recalculate mean and standard deviation.
        R
            for(i in 1:length(myFiles)){
                # read in table
                insertSizes <- read.table(file = myFiles[i])
                # output
                outFile <- gsub("NA128", "frag_size_dist/NA128", gsub(".txt", "_distribution_cor.png", basename(myFiles[i])))
                # sample name (for title)
                sampleName <- gsub(".sorted.insert_sizes.txt", "", basename(myFiles[i]))
                
                # Calculate cut-off for outliers
                quant <- quantile(insertSizes$V1, probs = c(0.25, 0.75))
                iqrange <- IQR(insertSizes$V1)
                outlierUp <- quant[2]+3*iqrange
                outlierLow <- quant[1]-3*iqrange
                # Correct dataset
                insertSizes <- filter(insertSizes, V1 > outlierLow & V1 < outlierUp)
                
                # mean and sd of insert sizes
                cat(c("\n", sampleName, "mean insert size:", "\t", mean(insertSizes$V1), "\n", sampleName, "SD of insert size:", "\t", sd(insertSizes$V1), sampleName, "\n", "summary statistics:", summary(insertSizes)))
                }

        q()

                #  NA12878_chr20 mean insert size: 	 -0.350630927561376 
                #  NA12878_chr20 SD of insert size: 	 326.876905874861 
                #  NA12878_chr20 summary statistics: 
                    # Min.   :-2165.0000   
                    # 1st Qu.: -319.0000   
                    # Median :  -50.0000   
                    # Mean   :   -0.3506   
                    # 3rd Qu.:  319.0000   
                    # Max.   : 2186.0000  

                #  NA12891_chr20 mean insert size: 	 -0.479521341123947 
                #  NA12891_chr20 SD of insert size: 	 307.816271697515 
                #  NA12891_chr20 summary statistics: 
                    # Min.   :-2052.0000   
                    # 1st Qu.: -300.0000   
                    # Median :  -76.0000   
                    # Mean   :   -0.4795   
                    # 3rd Qu.:  299.0000   
                    # Max.   : 2037.0000  

                #  NA12892_chr20 mean insert size: 	 -0.322018937203019 
                #  NA12892_chr20 SD of insert size: 	 312.294119759515
                #  NA12892_chr20 summary statistics: 
                    # Min.   :-1968.000   
                    # 1st Qu.: -303.000   
                    # Median :  -77.000   
                    # Mean   :   -0.322   
                    # 3rd Qu.:  302.000   
                    # Max.   : 2059.000  


    ## 9. Based on the file generated for sample NA12891 in 5, what would you say DELLY will use as a threshold for discordant mappings?
        ## hint: you will need to pipe the file into awk to print the column number corresponding to the MEDIAN_INSERT_SIZE and STANDARD_DEVIATION
        awk 'NR == 8 {print $1, $7}' frag_size_dist/NA12891_chr20.sorted.insert_size_metrics.txt
            ## MEDIAN_INSERT_SIZE   STANDARD_DEVIATION
            ## 300                  65.277149

            ## around 500 when calculated

# 4 SV detection
    ## Call germline SVs on each sample
    
    ## 10. The notation of the chromosome will cause an error. Change the first line in the reference file Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa from >20 to >chr20 using unixs sed command.
        sed -i '1 s/>20/>chr20/' refgenome/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa
        # index again
        samtools faidx refgenome/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa

        ## THIS IS VERY OFTEN THE PROBLEM! Some programs prefer with chr, some without

    ## 11. For each sample, type:
        for INPUT in data/NA128*_chr20.sorted.bam; do
            OUTPUT=$(echo $INPUT | sed "s/\.bam/\.bcf/" | sed "s/data\//svdetection\//")
            delly call -g refgenome/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa -o $OUTPUT \
                -x svdetection/human.hg38.excl.tsv \
                $INPUT
        done

    ## Look at the output you can use bcftools
        bcftools view svdetection/NA12878_chr20.sorted.bcf | less -S
        bcftools view svdetection/NA12891_chr20.sorted.bcf | less -S
        bcftools view svdetection/NA12892_chr20.sorted.bcf | less -S

    ## 12. generate a file with 
        touch delly_metrics.txt

    ## Pipe the following outputs of each bcf file into the file. Dont forget to use >> to append.
        for INPUT in svdetection/NA128*_chr20.sorted.bcf; do
            ## a) the file name
            basename $INPUT >> svdetection/delly_metrics.txt

            ## b) the number of SVs did you find in each sample
            bcftools view --no-header $INPUT | grep -v "^#" | wc -l >> svdetection/delly_metrics.txt

            ## c) the number of each type of SV (e.g., deletions, duplications, inversions) in each sample
            bcftools view $INPUT | grep -v "^#" | awk ' {print $5} ' | sort | uniq -c >> svdetection/delly_metrics.txt
        done

        ##DO THIS IN TABLE


# 4.1 Merging SV calls
    ## 13. Merge the SVs from multiple samples into a single, unified site list with DELLY:
        delly merge -m 500 -n 1000000 \
            -o svdetection/SVs_merged.bcf \
            -b 500 -r 0.5 svdetection/NA12878_chr20.sorted.bcf svdetection/NA12891_chr20.sorted.bcf svdetection/NA12892_chr20.sorted.bcf

    ## 14. Examine the output with bcftools. Pipe the following outputs of the merged bcf file into delly_metrics.txt.
        ## a) the file name
        basename svdetection/SVs_merged.bcf >> svdetection/delly_metrics.txt

        ## b) the number of SVs
        bcftools view --no-header svdetection/SVs_merged.bcf | grep -v "^#" | wc -l >> svdetection/delly_metrics.txt

        ## c) the number of each type of SV (e.g., deletions, duplications, inversions)
        bcftools view svdetection/SVs_merged.bcf | grep -v "^#" | awk ' {print $5} ' | sort | uniq -c >> svdetection/delly_metrics.txt


# 4.2 Re-genotyping
    ## 15. Genotype this merged SV site list across all samples.
    ## 16. Name the output .bcf files NA12878_geno.bcf, NA12891_geno.bcf, and NA12892_geno.bcf.
        for INPUT in data/NA128*_chr20.sorted.bam; do
            OUTPUT=$(echo $INPUT | sed "s/\_chr20\.sorted\.bam/\_geno\.bcf/" | sed "s/data\//svdetection\//")
            delly call -g refgenome/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa \
                -v svdetection/SVs_merged.bcf \
                -o $OUTPUT \
                -x svdetection/human.hg38.excl.tsv $INPUT
        done

    ## 17. Examine the output files with bcftools. You should see information on the genotype for each individual.
        bcftools view svdetection/NA12878_geno.bcf |  less -S
        bcftools view svdetection/NA12891_geno.bcf |  less -S
        bcftools view svdetection/NA12892_geno.bcf |  less -S
        

# 4.3 Merging the re-genotyped calls
    ## 18. Merge all re-genotyped samples to get a single .bcf output file using bcftools merge:
        bcftools merge -O b \
            -o svdetection/geno_merged.bcf \
            svdetection/NA12878_geno.bcf svdetection/NA12891_geno.bcf svdetection/NA12892_geno.bcf
            
        
    ## 19. Index the resulting .bcf file and create a .vcf file for visualization:
            bcftools index svdetection/geno_merged.bcf
            bcftools view svdetection/geno_merged.bcf > svdetection/geno_merged.vcf
            bgzip -c svdetection/geno_merged.vcf > svdetection/geno_merged.vcf.gz
            tabix -fp vcf svdetection/geno_merged.vcf.gz

# 5. Setting up IGV for SV visualization
    ## 20. Launch IGV and load the merged SV calls and the mapping files for the individual samples using File -> Load from File.
        svdetection/geno_merged.vcf.gz
        svdetection/geno_merged.vcf.gz.tbi
        data/NA12878_chr20.sorted.bam
        data/NA12891_chr20.sorted.bam
        data/NA12892_chr20.sorted.bam

    ## 21. Navigate to the following location to see a deletion: chr20:63,090,172-63,097,143.
    ## 22. Load the .bam files.
    ## You can try to configure IGV such that we can more clearly see the alignments that support the SV prediction(s).
    ## Colour the alignments by insert size and pair orientation.

# 5.1 Explore the SV
    ## 23. Is the variant at chr20:30,155,057-30,168,654 found in any member of the trio? If yes, how many reads support it?
        ## It's present in NA12891, although some reads still map to the junction. 13 reads support it.
                
                bcftools view svdetection/SVs_merged.bcf | grep -v '##' | head -1
                bcftools view svdetection/SVs_merged.bcf | grep "DEL00000031"
                # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
                # chr20   30160592        DEL00000031     G       <DEL>   504     PASS    IMPRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.9.1;END=30162407;PE=9;MAPQ=60;CT=3to5;CIPOS=-50,50;CIEND=-50,50

    ## 24. What are the genotypes for each member of the trio at the locus chr20:61,783,593-61,784,755 (e.g., hemizygous, homozygous)?
        ## NA12878  homozygous
        ## NA1291   homozygous
        ## NA1291   heterozygous

    ## 25. What about the variant at chr20:52,142,394-52,144,695?
        ## NA12878  homozygous for the variant
        ## NA1291   homozygous for the variant
        ## NA1291   heterozygous
