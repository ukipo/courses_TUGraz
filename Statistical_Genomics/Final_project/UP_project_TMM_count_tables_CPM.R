# Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("Mus.musculus")

library(dplyr)
library(limma)
# library(Glimma)
library(edgeR)
library(Mus.musculus)
library(biomaRt)


# Files with counts
files_raw <- Sys.glob(file.path("C:/Users/Uki/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/up0289/StatGen/project/counted_reads/raw_genome_tran/", 'ERR2588*.tsv'))

files_trim <- Sys.glob(file.path("C:/Users/Uki/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/up0289/StatGen/project/counted_reads/trimmed/", 'ERR2588*.tsv'))

files <- files_trim
files <- files_raw

# Import files
x <- readDGE(files, columns=c(1,2))  #CHECK WHICH COLUMNS
class(x)
dim(x)

# Organising sample information
# Rename samples
samplenames <- substring(colnames(x), (nchar(colnames(x)) - 32), (nchar(colnames(x)) - 23)) #CHECK WHICH FIRST, LAST
samplenames
colnames(x) <- samplenames
# Group samples by tissue
group <- as.factor(c("Liver", "Liver", "Brain", "Brain", "Brain", "Brain", "Liver", "Liver", "Brain", "Brain", "Liver", "Liver"))
x$samples$group <- group
# Group samples by tissue
timepoint <- as.factor(c("12.5", "12.5", "12.5", "12.5", "18.5", "18.5", "18.5", "18.5", "2wpb", "2wpb", "2wpb", "2wpb"))
x$samples$timepoint <- timepoint

# Data pre-processing

# Removing genes that are lowly expressed
  # some genes have no expression values across all samples
  table(rowSums(x$counts==0)==12)
    # TRIM: 36% of genes have zero counts across all nine samples
      # FALSE  TRUE 
      # 34929 19913
    # RAW: 36% of genes have zero counts across all nine samples
      # FALSE  TRUE 
      # 34864 19978
  
# Remove genes that did not reach at least 10 reads in both duplicates
  table(rowSums(x$counts>=10)==12)
  
# filter out genes with low or no expression
  # By default, the function keeps genes with about 10 read counts or more in a minimum number of samples, where the number of samples is chosen according to the minimum group sample size - 2. Uses CPM
  keep.exprs <- filterByExpr(x, group=group)
  x <- x[keep.exprs,, keep.lib.sizes=FALSE]
  dim(x)
  
# Transformations from the raw-scale
  # CPM
  cpm_x <- cpm(x)
  # LCPM
  lcpm_x <- cpm(x, log=TRUE)
  # RPKM
  # rpkm_x <- rpkm(x)
      # Error in rpkm.DGEList(x) : Gene lengths not found
  
# figure that shows the reduction of not expressed genes
  # Code to produce the figure is given below.
  #L and M
  L <- mean(x$samples$lib.size) * 1e-6
  M <- median(x$samples$lib.size) * 1e-6
  c(L, M)
# Figure
  x11()
pdf("Figure_filtered_data.pdf")
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm_x[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm_x[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm_x <- cpm(x, log=TRUE)
plot(density(lcpm_x[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm_x[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.off()

# Normalising gene expresison distributions - method of trimmed mean of M-values (TMM) (Robinson and Oshlack 2010)
  x <- calcNormFactors(x, method = "TMM")
  x$samples$norm.factors
  
# give a visual representation of the effects of normalisation
  x2 <- x
  x2$samples$norm.factors <- 1
  x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
  x2$counts[,2] <- x2$counts[,2]*5

# The figure below shows the expression distribution of samples for unnormalised and normalised data, where distributions are noticeably different pre-normalisation and are similar post-normalisation.
  x11()
pdf("Normalisation.pdf")
  par(mfrow=c(1,2))
  lcpm_x2 <- cpm(x2, log=TRUE)
  boxplot(lcpm_x2, las=2, col=col, main="")
  title(main="A. Example: Unnormalised data",ylab="Log-cpm")
  x2 <- calcNormFactors(x2)  
  x2$samples$norm.factors
  lcpm_x2 <- cpm(x2, log=TRUE)
  boxplot(lcpm_x2, las=2, col=col, main="")
  title(main="B. Example: Normalised data",ylab="Log-cpm")
dev.off()

# Export normalised gene counts
write.csv(cpm_x, file = "C:/Users/Uki/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/up0289/StatGen/project/normalised_readcounts/Mouse.trim.cpm.csv")
write.csv(lcpm_x, file = "C:/Users/Uki/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/up0289/StatGen/project/normalised_readcounts/Mouse.trim.lcpm.csv")

write.csv(cpm_x, file = "C:/Users/Uki/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/up0289/StatGen/project/normalised_readcounts/Mouse.raw.cpm.csv")
write.csv(lcpm_x, file = "C:/Users/Uki/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/up0289/StatGen/project/normalised_readcounts/Mouse.raw.lcpm.csv")
