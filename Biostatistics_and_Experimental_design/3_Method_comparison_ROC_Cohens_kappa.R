# 1. Method Comparison according to Bland-Altman [R]

# a.
# Read the dataset into R.
milk_fat <- read.csv("milkfat.csv")

# What can you say about the correlation of the measured data ("enzymatic_mehod" vs. "Gerber_method") with the aid of a scatter plot?
png(filename = "milk_scatter_plot.png")
plot(milk_fat$enzymatic_method, milk_fat$Gerber_method, type = "p", main = "Scatterplot of the enzymatic vs Gerber method", xlab="Enzymatic method", ylab="Gerber method", pch=19, col="orange", frame=FALSE)
abline(coef=c(0,1), col=2) # identity
abline(lm(milk_fat$enzymatic_method ~ milk_fat$Gerber_method), col="purple")
legend("topleft", legend=c("Identity", "Regression"), lwd=3, col=c(2, "purple"))


# Calculate the Pearson correlation coefficient to check your statement. What is the difference between correlation and agreement?
text(paste("Pearson Correlation: ", round(cor(milk_fat$enzymatic_method, milk_fat$Gerber_method, method="pearson"),3)), x=5, y=3, col="purple")
dev.off()

# b.
# Compare both fat content measurement methods, i.e. the measured data, using a Bland-Altman plot and interpret your results.
png(filename = "milk_bland_altman_plot.png")
x=(milk_fat$enzymatic_method + milk_fat$Gerber_method)/2
y=milk_fat$enzymatic_method - milk_fat$Gerber_method
plot(x, y, type = "p", main = "Bland-Altman plot comparrison of \nthe enzymatic vs Gerber method", xlab="(Enzymatic method + Gerber method)/2 [g/100 ml]", ylab="Enzymatic method - Gerber method [g/100 ml]", pch=19, col="blue4", frame=FALSE, xlim=c(0.5,7), ylim=c(-0.3, 0.3))
abline(h=mean(milk_fat$enzymatic_method - milk_fat$Gerber_method), col="maroon")
abline(h=mean(milk_fat$enzymatic_method - milk_fat$Gerber_method)+2*sd(milk_fat$enzymatic_method - milk_fat$Gerber_method), lty=2, col="maroon")
abline(h=mean(milk_fat$enzymatic_method - milk_fat$Gerber_method)-2*sd(milk_fat$enzymatic_method - milk_fat$Gerber_method), lty=2, col="maroon")
text(paste("\nmean"), x=max(x), y=mean(y), col="maroon")
text(paste("\nmean+2SD"), x=max(x), y=mean(y)+2*sd(y), col="maroon")
text(paste("\nmean-2SD"), x=max(x), y=mean(y)-2*sd(y), col="maroon")
dev.off()


# Sum up your commands for the construction of the Bland-Altman plot in a function and apply it. (Do not use the available R-package for the construction of the Bland-Altman plot. Write your own code and wrap it in a user-defined function.)
bland.altman <- function(methodA, methodB){
  x <- (methodA + methodB)/2
  y <- (methodA - methodB)
  nameA <- sub(".*\\$", "", deparse(substitute(methodA)))
  nameB <- sub(".*\\$", "", deparse(substitute(methodB)))
  plot(x, y, type = "p", main = paste("Bland-Altman plot comparrison of\n", nameA, " and ", nameB), xlab=paste("(", nameA, " + ", nameB, ")/2   [g/100 ml]"), ylab=paste(nameA, " - ", nameB, "   [g/100 ml]"), pch=19, col="blue4", frame=FALSE)
  abline(h=mean(y), col="maroon")
  abline(h=mean(y)+2*sd(y), lty=2, col="maroon")
  abline(h=mean(y)-2*sd(y), lty=2, col="maroon")
  text(paste("\nmean"), x=max(x), y=mean(y), col="maroon")
  text(paste("\nmean+2SD          "), x=max(x), y=mean(y)+2*sd(y), col="maroon")
  text(paste("\nmean-2SD          "), x=max(x), y=mean(y)-2*sd(y), col="maroon")
}

png(filename = "milk_function_ba_plot.png")
bland.altman(milk_fat$enzymatic_method, milk_fat$Gerber_method)
dev.off()

# c.
# Bonus task: Generate the Bland-Altman plot for the fat content data with the R-package BlandAltmanLeh.
png(filename = "milk_BlandAltmanLeh_plot.png")
library(BlandAltmanLeh)
bland.altman.plot(milk_fat$enzymatic_method, milk_fat$Gerber_method, main = paste("Bland-Altman plot comparrison of \nthe enzymatic and Gerber methods"), xlab="enzymatic + Gerber method)/2   [g/100 ml]", ylab="enzymatic - Gerber method   [g/100 ml]", pch=19, col="blue4", frame=FALSE)
dev.off()

# 2. ROC curve and diagnostic tests

# a.
# The test is positive if the blood sugar value is greater than the chosen cutoff, regardless whether the patient is diseased or not. Choose a cutoff of 5 and two other suitable cutoffs and calculate SN, SP, FPR, NPV, prevalence p and accuracy.
blood_sugar <- read.csv("blood_sugar.csv")
blood_sugar$lim <- c(5,6,7,7.1)

roc_table <- data.frame(matrix(ncol=12, nrow=3))
colnames(roc_table) <- c("cutoff", "TP", "FP", "FN", "TN", "SN", "SP", "FPR", "PPV", "NPV", "prevalence", "accuracy")
roc_table$cutoff=c(5, 6, 7)


for(i in 1:nrow(roc_table)){
  co <- roc_table$cutoff[i]
  roc_table$FN[i] <- sum(blood_sugar$sick[c(which(blood_sugar$lim <= co))])
  roc_table$TN[i] <- sum(blood_sugar$healthy[c(which(blood_sugar$lim <= co))])
  roc_table$TP[i] <- sum(blood_sugar$sick[c(which(blood_sugar$lim > co))])
  roc_table$FP[i] <- sum(blood_sugar$healthy[c(which(blood_sugar$lim > co))])
  roc_table$SN[i] <- roc_table$TP[i]/(roc_table$TP[i]+roc_table$FN[i])
  roc_table$SP[i] <- roc_table$TN[i]/(roc_table$FP[i]+roc_table$TN[i])
  roc_table$FPR[i] <- 1-roc_table$SP[i]
  roc_table$PPV[i] <- roc_table$TP[i]/(roc_table$FP[i]+roc_table$TP[i])
  roc_table$NPV[i] <- roc_table$TN[i]/(roc_table$TN[i]+roc_table$FN[i])
  n <- sum(roc_table$TP[i], roc_table$FP[i], roc_table$FN[i], roc_table$TN[i])
  roc_table$prevalence[i] <- (roc_table$TP[i]+roc_table$FN[i])/n
  roc_table$accuracy[i] <- (roc_table$TP[i]+roc_table$TN[i])/n
}

# b.
# Draw the ROC Curve in R and determine the most suitable cutoff for the diagnosis. Justify your decision!
png(filename = "blood_sugar_ROC.png")
plot(c(1, roc_table$FPR, 0), c(1, roc_table$SN, 0), type = "b", pch = 19, col = "red", xlab = "FPR", ylab = "SN", xlim = c(0,1), ylim = c(0,1))
abline(coef=c(0,1))
text(roc_table$FPR+0.02, roc_table$SN-0.02, labels = roc_table$cutoff)
dev.off()

# c.
# Assume you do a ROC curve, but instead of above the reference line, your curve appears below it. What is a possible explanation?


# 3. Cohen's kappa [R]

# a.
# Calculate the kappa statistics (unweighted, linear and square weighted) for the values given in Table 2 without using a predefined R package. Interpret your results. (Hint: You will need a double loop.)

# data
paptest <- data.frame(PA=c(rep("1", 45), rep("2", 35), rep("3", 42), rep("4", 23), rep("5", 21)),
                      PB=c(rep("1", 31), rep("2", 11), rep("3", 0), rep("4", 3), rep("5", 0), rep("1", 4), rep("2", 18), rep("3", 1), rep("4", 10), rep("5", 2), rep("1", 10), rep("2", 4), rep("3", 24), rep("4", 4), rep("5", 0), rep("1", 0), rep("2", 1), rep("3", 9), rep("4", 12), rep("5", 1), rep("1", 4), rep("2", 0), rep("3", 2), rep("4", 0), rep("5", 15)))

paptest <- table(paptest)

# Kappa
cohens_kappa <- function(x){
  # Observed and expected frequencies
  obs_freq <- x/sum(x)
  exp_freq <- colSums(x) %o% rowSums(x)/sum(x)^2
  # Unweighted Kappa
  pobs <- sum(diag(obs_freq))
  pexp <- sum(diag(exp_freq))
  uw_kappa <- (pobs-pexp)/(1-pexp)
  # Weight matrix
  weigth_linear <- matrix(nrow=nrow(x), ncol=ncol(x))
  weigth_square <- matrix(nrow=nrow(x), ncol=ncol(x))
  for (i in 1:nrow(weigth_linear)){
    for (j in 1:ncol(weigth_linear)){
      weigth_linear[i, j] <- 1-(abs(i-j)/(nrow(x)-1))
      weigth_square[i, j] <- 1-((abs(i-j))^2/(nrow(x)-1)^2)
    }
  }
  # Linear weighted Kappa
  lw_pobs <- sum(weigth_linear*obs_freq)
  lw_pexp <- sum(weigth_linear*exp_freq)
  lw_kappa <- (lw_pobs-lw_pexp)/(1-lw_pexp)
  # Square weighted Kappa
  sq_pobs <- sum(weigth_square*obs_freq)
  sq_pexp <- sum(weigth_square*exp_freq)
  sq_kappa <- (sq_pobs-sq_pexp)/(1-sq_pexp)
  # Return
  final <- data.frame("Unweighted_Kappa"=uw_kappa,
                      "Linear_weighted_Kappa"=lw_kappa,
                      "Square_weighted_Kappa"=sq_kappa)
  return(final)
}

cohens_kappa(paptest)


# b.
# Use the R package vcd to calculate the kappa statistics and compare your results.
library(vcd)
Kappa(paptest)
Kappa(paptest, weights = "Fleiss-Cohen")

# c.
# Considering the example from above, which of the kappa coefficients (unweighted, linear or square weighted) would you choose to describe the agreement and why?

