# Libraries

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("vsn")

library(DESeq2)
library(vsn)
library(FactoMineR)


  # Mouse data
  trimCts <- Sys.glob(file.path("~/StatGen/project/normalised_readcounts/Mouse.trim.cpm.csv"))

  trimAnno <- Sys.glob(file.path("~/StatGen/project/normalised_readcounts/trim_sample_annotation.csv"))
  
  cts_trim <- as.matrix(read.csv(trimCts,row.names="gene_id"))
  coldata_trim <- read.csv(trimAnno, row.names=1)
  coldata_trim <- coldata_trim[,c("condition","type")]
  coldata_trim$condition <- factor(coldata_trim$condition)
  coldata_trim$type <- factor(coldata_trim$type)


# IMPORT DATA FROM HTSEQ-COUNT ANALYSIS 
  
# Set directory of input files
  directory <- "~/StatGen/project/counted_reads/trimmed/"
  
# Specify which files to read with list.files
  sampleFiles <- grep("ERR2588",list.files(directory), value=TRUE)
  sampleNames <-  sub("*(.trim).sorted.readcounts.tsv","",sampleFiles)

  sampleCondition <- coldata_trim$condition
  sampleTissue <- factor(c("Liver", "Liver", "Brain", "Brain", "Brain", "Brain", "Liver", "Liver", "Brain", "Brain", "Liver", "Liver"), levels = c("Brain", "Liver"))
  # Timepoint
  sampleTimepoint <- c("12.5","12.5", "12.5", "12.5", "18.5", "18.5", "18.5", "18.5", "32", "32", "32", "32")
  
  sampleTable <- data.frame(sampleName = sampleNames,
                            fileName = sampleFiles,
                            condition = sampleCondition,
                            tissue = sampleTissue,
                            timepoint = sampleTimepoint)
  sampleTable$condition <- factor(sampleTable$condition)

# Build DESeqDataSet
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                         directory = directory,
                                         design= ~ c(condition, tissue, timepoint))
  ddsHTSeq
  
# Factor levels
  ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = c("Brain_12_5", "Liver_12_5", "Brain_18_5", "Liver_18_5", "Brain_2wpb", "Liver_2wpb"))
  
# Removing genes that are lowly expressed
  # some genes have no expression values across all samples
  table(rowSums(counts(ddsHTSeq)==0)<12)

  # discard genes that have no reads in any of the samples
  keep <- rowSums(counts(ddsHTSeq)==0)<12
  ddsHTSeq <- ddsHTSeq[keep,]
  
        # # Differential expression analysis
        #   dds <- DESeq(ddsHTSeq)
        #   res <- results(dds)  
        #   res <- results(dds, contrast=c("condition","Brain_12_5","Brain_2wpb"))
  
  
# VST NORMALIZATION
  
  # produce transformed data on the log2 scale which has been normalized with respect to library size or other normalization factors
  
  # Blind dispersion estimation
    # use defaulte blind = TRUE, so that the comparison is unbiased (differences in counts are not assumed to be explainable by the experimental design)
  
# Extracting transformed values - variance stabilizing transformation is used
  vsd <- vst(ddsHTSeq)
  head(assay(vsd), 3)
  head(assay(ddsHTSeq), 3)
  
# Visualize the stabilized variance - sd should be roughly constant across the whole dynamic range
  ntd <- normTransform(ddsHTSeq)
  library("vsn")
  x11()
  meanSdPlot(assay(ntd))
  meanSdPlot(assay(vsd))
  
  
# PCA ANALYSIS (FactoMineR)
  
# Define conditions
  # Tissue
  vsd@colData@listData[["tissue"]] <- factor(c("Liver", "Liver", "Brain", "Brain", "Brain", "Brain", "Liver", "Liver", "Brain", "Brain", "Liver", "Liver"), levels = c("Brain", "Liver"))
  # Timepoint
  vsd@colData@listData[["timepoint"]] <- factor(c("12_5","12_5", "12_5", "12_5", "18_5", "18_5", "18_5", "18_5", "2wpb", "2wpb", "2wpb", "2wpb"), levels = c("12_5", "18_5", "2wpb"))
  
  
# The dist function expects the different samples to be rows of its argument, and different dimensions - genes - to be columns
  tvsd <- t(assay(vsd))
  
# Side quest: Look at the heat map
  library(pheatmap)
  sampleDists <- dist(tvsd, method = "euclidean")
  x11()
  pheatmap(sampleDists, labels_row =sampleCondition)
  
# PCA
  pca.data <- PCA(tvsd, graph = FALSE)
  
  # plot
  library(ggplot2)
  x11()
  plot(pca.data, choix="ind", col=)
  plotPCA(vsd, intgroup="condition")
  plotPCA(vsd, intgroup="tissue")
  plotPCA(vsd, intgroup="timepoint")
  plotPCA(vsd, intgroup=c("tissue", "timepoint"))
  
  plot.PCA(pca.data, choix="ind", graph.type = "ggplot")
  
  # sampleDistMatrix <- as.matrix(sampleDists)
  # rownames(sampleDistMatrix) <- paste(vsd$tissue, vsd$timepoint, sep="_")
  # colnames(sampleDistMatrix) <- NULL
  # colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  # pheatmap(sampleDistMatrix,
  #          clustering_distance_rows=sampleDists,
  #          clustering_distance_cols=sampleDists,
  #          col=colors)
  
  plotPCA(vsd, intgroup=c("tissue", "timepoint"))
  
  pcaData <- plotPCA(vsd, intgroup=c("tissue", "timepoint"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p <- ggplot(pcaData, aes(PC1, PC2, color=tissue, size=timepoint, col)) +
    geom_point() +
    scale_color_manual(values=c("dodgerblue3", "green3")) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
  
  p + theme_light() + scale_size_discrete(range = c(3.5, 6.5, 10))
  
  
# DEA with DESeq
  sampleTable$time <- rep(c(12.5, 18.5, 32), each=4)
  
  dds <- DESeq(ddsHTSeq, test = "Wald", fitType = "parametric")
  ddsHTSeq@assays@data@listData <- as.vector(ddsHTSeq@assays@data@listData)
  dds <- DESeq(ddsHTSeq)
  
  res <- results(ddsHTSeq)
