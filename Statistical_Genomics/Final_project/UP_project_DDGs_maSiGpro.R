# Libraries
  library(maSigPro)
  library(dplyr)
  library(tidyr)
  library(ggplot2)

# Data - count tables from EdgeR (in CPM)
  data <- read.csv("Mouse.trim.cpm.csv")
  # In each organ
  data_brain <- dplyr::select(data, "ERR2588422", "ERR2588424", "ERR2588545", "ERR2588547", "ERR2588592", "ERR2588593")
  data_liver <- dplyr::select(data, "ERR2588414", "ERR2588416", "ERR2588561", "ERR2588563", "ERR2588607", "ERR2588608")
  
# Experimental design
  edesign <- read.csv("edesign_mouse_trim.txt")
  edesign_brain <- read.csv("edesign_mouse_trim_brain.csv")
  edesign_liver <- read.csv("edesign_mouse_trim_liver.csv")
  
# log-transformed time (measured in days psot-conception)
  edesign$Time <- log(edesign$Time)
  edesign_brain$Time <- log(edesign_brain$Time)
  edesign_liver$Time <- log(edesign_liver$Time)

# Regression matrix - degree = 3
  design <- make.design.matrix(edesign, degree = 3)
  design_brain <- make.design.matrix(edesign_brain, degree = 3)
  design_liver <- make.design.matrix(edesign_liver, degree = 3)
  
  ss.design <- make.design.matrix(ss.edesign)
  
# Compute a regression fit for each gene
  fit <- p.vector(data=data, design=design, counts = TRUE)
  fit_brain <- p.vector(data=data_brain, design=design_brain, counts = TRUE)
  fit_liver <- p.vector(data=data_liver, design=design_liver, counts = TRUE)
  
  ss.fit <- p.vector(ss.DATA, ss.design, counts = TRUE)
  
# Find significant variables for each gene
  tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
  tstep_brain <- T.fit(fit_brain, step.method = "backward", alfa = 0.05)
  tstep_liver <- T.fit(fit_liver, step.method = "backward", alfa = 0.05)
  
  ss.tstep <- T.fit(ss.fit, step.method = "backward", alfa = 0.05)
  
# Generate lists of significant genes
  sigs <- get.siggenes(tstep = tstep, rsq = 0.3, vars = "each")
  sigs_all <- get.siggenes(tstep = tstep, rsq = 0.3, vars = "all")
  head(names(sigs))
  head(names(sigs_all))
  head(sigs$summary)
  head(sigs_all$summary)
  
  ss.sigs <- get.siggenes(ss.tstep, rsq = 0.3, vars = "each")
  
  # each
  sigs_brain <- get.siggenes(tstep = tstep_brain, rsq = 0.3, vars = "each")
  head(sigs_brain$summary)
  
  sigs_liver <- get.siggenes(tstep = tstep_liver, rsq = 0.3, vars = "each")
  head(sigs_liver$summary)
  
  # all
  sigs_brain_all <- get.siggenes(tstep = tstep_brain, rsq = 0.3, vars = "all")
  head(sigs_brain$summary)
  
  sigs_liver_all <- get.siggenes(tstep = tstep_liver, rsq = 0.3, vars = "all")
  head(sigs_liver$summary)
  
# Graphic display
  suma2Venn(sigs_brain$summary[, c(1:4)])
  suma2Venn(sigs_liver$summary[, c(1:4)])
  
# Number of differentially expressed genes
  sigs_brain$sig.genes$independ$g
  sigs_brain$sig.genes$Time$g
  sigs_brain$sig.genes$Time2$g
  sigs_brain$sig.genes$Time3$g
  
  sigs_liver$sig.genes$independ$g
  sigs_liver$sig.genes$Time$g
  sigs_liver$sig.genes$Time2$g
  sigs_liver$sig.genes$Time3$g
  
# Count genes
  
  expr_genes <- data.frame(matrix(ncol=10, nrow=length(data_brain$ERR2588422)), row.names = rownames(data_brain))
  colnames(expr_genes) <- c("Brain_independent", "Brain_Time", "Brain_Time2", "Brain_Time3", "Liver_independent", "Liver_Time", "Liver_Time2", "Liver_Time3", "Brain", "Liver")
  
  
  for(i in 1:length(sigs_brain$summary$independ)){
    n <- which(rownames(expr_genes)==sigs_brain$summary$independ[i])
    expr_genes$Brain_independent[n]=1
    expr_genes$Brain[n]=1
    m <- which(rownames(expr_genes)==sigs_brain$summary$Time[i])
    expr_genes$Brain_Time[m]=1
    expr_genes$Brain[m]=1
    l <- which(rownames(expr_genes)==sigs_brain$summary$Time2[i])
    expr_genes$Brain_Time2[l]=1
    expr_genes$Brain[l]=1
    o <- which(rownames(expr_genes)==sigs_brain$summary$Time3[i])
    expr_genes$Brain_Time3[o]=1
    expr_genes$Brain[o]=1
    rm(n)
    rm(m)
    rm(l)
    rm(o)
  }
  
  for(i in 1:length(sigs_liver$summary$independ)){
    n <- which(rownames(expr_genes)==sigs_liver$summary$independ[i])
    expr_genes$Liver_independent[n]=1
    expr_genes$Liver[n]=1
    m <- which(rownames(expr_genes)==sigs_liver$summary$Time[i])
    expr_genes$Liver_Time[m]=1
    expr_genes$Liver[m]=1
    l <- which(rownames(expr_genes)==sigs_liver$summary$Time2[i])
    expr_genes$Liver_Time2[l]=1
    expr_genes$Liver[l]=1
    o <- which(rownames(expr_genes)==sigs_liver$summary$Time3[i])
    expr_genes$Liver_Time3[o]=1
    expr_genes$Liver[o]=1
    rm(n)
    rm(m)
    rm(l)
    rm(o)
  }

  expr_genes[is.na(expr_genes)] <- 0
  expr_genes <- expr_genes[1:16072,]
  
  for (i in 1:ncol(expr_genes)) {
    expr_genes[16073,i] <- sum(expr_genes[1:16072,i])
  }
  
  write.table(expr_genes, file = "DDGs.txt")
  
# Plots
  # subset
  no_ddg <- as.data.frame(t(expr_genes[16073,])); no_ddg
  no_ddg <- dplyr::slice(no_ddg, c(2:4, 6:8)); no_ddg
  no_ddg$time <- c(12.5, 18.5, 32, 12.5, 18.5, 32); no_ddg
  no_ddg$tissue <- rep(c("Brain", "Liver"), each=3); no_ddg
  
  
  ggplot(no_ddg, aes(time, sum, colour = tissue)) +
    facet_grid(cols = vars(tissue)) +
    geom_line()
  
  
  data_mouse <- read.table("Mouse_CPM.txt", header=TRUE, sep=" ")
  data_mouse <- dplyr::select(data_mouse, Brain.e12.5.1, Brain.e12.5.2, Brain.e12.5.3, Brain.e12.5.4, Brain.e18.5.1, Brain.e18.5.2, Brain.e18.5.3, Brain.e18.5.4, Brain.P14.1, Brain.P14.2, Brain.P14.3, Brain.P14.4)
  edesign_mouse <- read.csv("edesign_mouse_cpm.csv")
  edesign_mouse$Time <- log(edesign_mouse$Time)
  design_mouse <- make.design.matrix(edesign_mouse, degree = 3)
  fit_mouse <- p.vector(data=data_mouse, design=design_mouse, counts = TRUE)
  tstep_mouse <- T.fit(fit_mouse, step.method = "backward", alfa = 0.05)
  sigs_mouse <- get.siggenes(tstep = tstep_mouse, rsq = 0.3, vars = "each")
  
  expr_mouse <- data.frame(matrix(ncol=4, nrow=length(data_mouse$Brain.e12.5.1)), row.names = rownames(data_mouse))
  colnames(expr_mouse) <- c("Brain_independent", "Brain_Time", "Brain_Time2", "Brain_Time3")
  
  
  for(i in 1:length(sigs_mouse$summary$independ)){
    n <- which(rownames(expr_mouse)==sigs_mouse$summary$independ[i])
    expr_mouse$Brain_independent[n]=1
    m <- which(rownames(expr_mouse)==sigs_mouse$summary$Time[i])
    expr_mouse$Brain_Time[m]=1
    l <- which(rownames(expr_mouse)==sigs_mouse$summary$Time2[i])
    expr_mouse$Brain_Time2[l]=1
    o <- which(rownames(expr_mouse)==sigs_mouse$summary$Time3[i])
    expr_mouse$Brain_Time3[o]=1
    rm(n)
    rm(m)
    rm(l)
    rm(o)
  }
  
  expr_mouse[is.na(expr_mouse)] <- 0
  
  for (i in 1:ncol(expr_mouse)) {
    expr_mouse[6437,i] <- sum(expr_mouse[1:6436,i])
  }
  