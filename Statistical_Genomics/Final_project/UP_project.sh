#!/bin/bash

# Foreword
# This project aims to partially replicate the results of the study published in the following paper:
#   Cardoso-Moreira, M., Halbert, J., Valloton, D. et al. Gene expression across mammalian organ development. Nature 571, 505–509 (2019). https://doi.org/10.1038/s41586-019-1338-505–509
# The data was obtained from the paper's repository and the described methods were adapted from said paper.


# 0 Data

    # Data files
    mkdir -p data
    cd data
        ## 1) connect to the ensembl server:
        ftp ftp.sra.ebi.ac.uk

        ## 2) navigate to the correct directory with the cd command
        cd vol1/fastq/ERR258/008/ERR2588568/

        ## 3) download the assembly with the get command
        declare -a StringArray=("ERR2588414/ERR2588414.fastq.gz" "ERR2588545/ERR2588545.fastq.gz" "ERR2588592/ERR2588592.fastq.gz" "ERR2588416/ERR2588416.fastq.gz" "ERR2588547/ERR2588547.fastq.gz" "ERR2588593/ERR2588593.fastq.gz" "ERR2588422/ERR2588422.fastq.gz" "ERR2588561/ERR2588561.fastq.gz" "ERR2588607/ERR2588607.fastq.gz" "ERR2588424/ERR2588424.fastq.gz" "ERR2588563/ERR2588563.fastq.gz" "ERR2588608/ERR2588608.fastq.gz")

        for i in SAMPLENAMES; do
            get $i
        done

        quit

    ls data/


    # Mouse Mus musculus data
    # 3 time points: E12.5, E18.5, P14 (2wpb)
    # 2 tissues: Brain, Liver
    # 2 replicates

    # Reference genome
    # Mus musculus GRCm38 (mm10)
    mkdir -p refgenome
        ## 1) connect to the ensembl server:
        ftp ftp.ensembl.org

        ## 2) navigate to the correct directory with the cd command
        ## Since the repository was last update on 18.02.2019, I used the release 95, which was the most recent version at the time (09.01.2019)
        cd pub/release-95/fasta/mus_musculus/dna/

        ## 3) download the assembly with the get command
        get Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz

        ## 4) navigate to the gtf directory
        cd ../../../gtf/mus_musculus/

        ## 5) download the position of exons
        Mus_musculus.GRCm38.95.gtf.gz

        ## 6) quit
        quit

    # Protocol description
    # RNA was extracted using the RNeasy protocol from QIAGEN. RNeasy Micro columns were used to extract RNA from small (< 5 mg) or fibrous samples and RNeasy Mini columns were used to extract RNA from larger samples. The tissues were homogenized in RLT buffer supplemented with 40 mM dithiothreitol (DTT) or QIAzol.	The organs were dissected from 14 developmental stages (starting at e10.5 until P63). Until e12.5 the dissected â€˜brainâ€™ samples correspond to the whole brain, and from e13.5 onwards they correspond to the cerebral hemispheres (without the olfactory bulbs). Early â€˜hindbrainâ€™ samples correspond to the prepontine hindbrain-enriched brain region (until e15.5) and from then onwards only to the cerebellum.	The RNA-seq libraries were sequenced on the HiSeq 2500 platform and were multiplexed in sets of 6 or 8.	The RNA-seq libraries were created using the TruSeq Stranded mRNA LT Sample Prep Kit (Illumina).	We normalized the count data using the method TMM as implemented in the package EdgeR (3.14.0). EdgeR was also used to generate the expression tables. Expression levels were calculated as cpm (counts per million) or in RPKM (reads per kilobase of exon model per million mapped reads)	We mapped the reads from each library against the reference GRCm38 (mm10) using GSNAP (22-10-2014). The alignments were guided by the known gene annotations (Ensembl E71) and the discovery of novel splice sites was enabled

# 1 Quality control with FastQC

    # Make directory to store analysis
    mkdir -p fastqc

    ## Run FastQC for each .fastq file
    for INPUT in data/ERR2588*.fastq.gz; do
        fastqc -o fastqc --nogroup $INPUT
    done

# 2.1 Trim the reads

    # Trim
    for INPUT in data/ERR2588*.fastq.gz; do
        OUTPUT=$(echo $INPUT | sed "s/\.fastq\.gz/\.trim\.fastq\.gz/" | sed "s/data\//trimmed\//")
        java -jar ~/progs/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 2 -phred33 $INPUT $OUTPUT ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:36 AVGQUAL:20
    done

    ## Run FastQC for each .trim.fastq file
    for INPUT in data/ERR2588*.trim.fastq.gz; do
        fastqc -o fastqc --nogroup $INPUT
    done

# 2.2 Prepare data

    # make exon, splicesite file
    hisat2_extract_splice_sites.py refgenome/Mus_musculus.GRCm38.95.gtf > refgenome/Mus_musculus.GRCm38.95.ss
    hisat2_extract_exons.py refgenome/Mus_musculus.GRCm38.95.gtf > refgenome/Mus_musculus.GRCm38.95.exon

    # Build HGFM index with transcripts
    hisat2-build -p 16 --exon refgenome/Mus_musculus.GRCm38.95.exon --ss refgenome/Mus_musculus.GRCm38.95.ss refgenome/Mus_musculus.GRCm38.dna_sm.toplevel.fa refgenome/genome_tran

# 2.3 Mapping and counting with Hissat2

    # Make directory to store analysis
    mkdir -p mapped_reads
    
    # Align trimmed reads
    for INPUT in trimmed/ERR2588*.trim.fastq.gz; do
        echo "aligning sample $INPUT"
        OUTPUT=$(echo $INPUT | sed "s/\.fastq\.gz/\.bam/" | sed "s/data\//mapped_reads\//")
        hisat2 -q -x refgenome/genome_tran -U $INPUT | samtools view -Sbo $OUTPUT -
        echo "done with sample $INPUT"
    done

    # Align raw reads
    for INPUT in data/ERR2588*.fastq.gz; do
        echo "aligning sample $INPUT"
        OUTPUT=$(echo $INPUT | sed "s/\.fastq\.gz/\.bam/" | sed "s/data\//mapped_reads\//")
        hisat2 -q -x refgenome/genome_tran -U $INPUT | samtools view -Sbo $OUTPUT -
        echo "done with sample $INPUT"
    done

    # Sort and index files with samtools
    for INPUT in mapped_reads/ERR2588*.bam; do
        echo "sorting and indexing file $INPUT"
        OUTPUT=$(echo $INPUT | sed "s/\.bam/\.sorted\.bam/")
        samtools sort -O 'bam' -T tmp_ -o $OUTPUT $INPUT
        samtools index -b $OUTPUT
    done


# 3 PCA on normalized (VST) counts (DESeq2 and R's prcomp() function)
    # Generate read counts for the set of protein-coding genes (HTSeq v.0.6.1)
            # # https://htseq.readthedocs.io/en/master/htseqcount.html?highlight=count
            # # htseq-count: counting reads within features
            # htseq-count [options] <alignment_files> <gtf_file> -s reverse
            # python -m HTSeq.scripts.count instead of htseq-count if "command not found"

        # Make directory
        mkdir -p counted_reads

        # Count reads per feature
        for INPUT in mapped_reads/trimmed/ERR2588*.sorted.bam; do
            echo "counting reads per feature in file $INPUT"
            OUTPUT=$(echo $INPUT | sed "s/\.bam/\.readcounts\.tsv/" | sed "s/mapped_reads\//counted_reads\//")
            htseq-count -f bam -r pos -s yes -c $OUTPUT $INPUT refgenome/Mus_musculus.GRCm38.95.gtf
        done

        for INPUT in mapped_reads/raw_genome_tran/ERR2588*.sorted.bam; do
            echo "counting reads per feature in file $INPUT"
            OUTPUT=$(echo $INPUT | sed "s/\.bam/\.readcounts\.tsv/" | sed "s/mapped_reads\//counted_reads\//")
            htseq-count -f bam -r pos -s yes -c $OUTPUT $INPUT refgenome/Mus_musculus.GRCm38.95.gtf
        done

        htseq-count -f bam -r pos -s yes -c counted_reads/Mouse.readcounts.tsv mapped_reads/trimmed/ERR2558*sorted.bam refgenome/Mus_musculus.GRCm38.95.gtf
    
    # Variance stabilizing transformation (VST) (DESeq2 v.1.12.4)
            # DESeq2 vignette
            # Data transformation and visualization
            # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

        # R SCRIPT
        # UP_project_VST_PCA.R chapter VST NORMALIZATION


        
    # Evaluate seq libraries with unsupervised hierarchical clustering (hclust) and PCA (FactoMineR 1.34) - input is read counts after VST
            # FactoMineR PCA description
            # http://factominer.free.fr/factomethods/principal-components-analysis.html

            # Principal Component Methods in R: Practical Guide
            # http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/117-hcpc-hierarchical-clustering-on-principal-components-essentials/

            # FactoMineR: An R Package for Multivariate Analysis
            # with included R script
            # https://www.jstatsoft.org/article/view/v025i01

        # R SCRIPT
        # UP_project_VST_PCA.R chapter PCA ANALYSIS (FactoMineR)

    # Generate expression tables used in the study in counts per million (CPM) or reads per kilobase of exon model per million mapped reads (RPKM) (EdgeR)
            # https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
            # 5 Data pre-processing
            # 5.1 Transformation from the raw sc

        # R SCRIPT
        # UP_project_TMM_count_tables_CPM.R


    # Normalize count data using the method TMM as implemented in EdgeR (v.3.14.0)
            # https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
            # Workflow: Normalising gene expression distributions

        # R SCRIPT
        # UP_project_TMM_count_tables_CPM.R


    # Correlation among replicates (Spearman's p) > 0.90
        # all samples had above, the others were excluded from avaliable samples


    # Identify DDGs using maSigPro
            # maSigPro Users Guide
            # https://www.bioconductor.org/packages/release/bioc/vignettes/maSigPro/inst/doc/maSigProUsersGuide.pdf
            # input=count tables from EdgeR in CPM, exclude genes that did not reach a minimum of 10 reads in (at least 3 libraries) both duplicates
            # regression matrix degree = 3

        # R SCRIPT
        # UP_project_DDGs_maSiGpro.R
