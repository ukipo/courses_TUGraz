# 4th Exercise WS 2021/22

# 1. Transformations [R]
# In a study about traffic safety the blood alcohol content (BAC) of 413 drunk drivers was measured. The legal limit for driving is at 0.05%, which equals 50 mg/dL BAC [1]. The according data file "drunk_driving.csv" can be downloaded on the course homepage.
# For further statistical analysis it is useful to have normally distributed data. For this reason, several transformation techniques can be used. This task should be solved in R.
  
  # import data
  drunk_driving <- read.csv("drunk_driving.csv", row.names = 1)

# [0.5P] a.
# Generate a qq-plot and a histogram to check if the given dataset follows a normal distribution and interpret your plots.
  
  # qq plot
  png(filename = "drunk_driving_qqplot.png")
    qqnorm(drunk_driving$BAC, main="Normal Q-Q Plot of blood alcohol content of drunk drivers", col = "violetred4")
    qqline(drunk_driving$BAC)
  dev.off()
  
  # histogram
  png(filename = "drunk_driving_hist.png")
    hist(drunk_driving$BAC, main="Blood alcohol content of drunk drivers", xlab="Blood alcohol content [mg/dL]", col = "slateblue1")
  dev.off()

# [2P] b.
# Transform the data by using the logarithm dualis, the square root transformation and the reciprocal transformation. Generate histograms and qq-plots of the transformed data and display them in two figures (one for histograms, one for qq-plots). Looking at the plots, which transformation works best? Justify your answer.
  
  # logarithm dualis transform
  dd_log <- log2(drunk_driving$BAC)
  
  # square root transform
  dd_sq <- sqrt(drunk_driving$BAC)
  
  # reciprocal transform
  dd_rt <- 1/drunk_driving$BAC
  
  # plot histograms
  png(filename = "drunk_driving_transform_hist.png") # open: write .png
    par(mfrow = c(2,2)) #parameters - 2 rows 2 columns for plots
    # untransformed
    hist(drunk_driving$BAC, main="Blood alcohol content of drunk drivers\nUntransformed data", xlab="Blood alcohol content [mg/dL]", col = "slateblue1")
    
    # logtransformed
    hist(dd_log, main="Blood alcohol content of drunk drivers\nLogarithm dualis transformation", xlab="Log2 of blood alcohol content [mg/dL]", col = "skyblue")
    
    # square root transformed
    hist(dd_sq, main="Blood alcohol content of drunk drivers\nSquare root transformation", xlab="Square root of blood alcohol content [mg/dL]", col = "skyblue")
    
    # reciprocal transformed
    hist(dd_rt, main="Blood alcohol content of drunk drivers\nReciprocal transformation", xlab="Reciprocal value of blood alcohol content [mg/dL]", col = "skyblue")
      
    par(mfrow = c(1,1)) # restore parameters
  dev.off() # close: write png
  
  # plot qq-plots
  png(filename = "drunk_driving_transform_qqplot.png")
    par(mfrow = c(2,2)) #parameters - 2 rows 2 columns for plots
    # untransformed
    qqnorm(drunk_driving$BAC, main="Blood alcohol content of drunk drivers\nqq plot - Untransformed data", col = "violetred4")
    qqline(drunk_driving$BAC)
    
    # logtransformed
    qqnorm(dd_log, main="Blood alcohol content of drunk drivers\nqq plot - Logarithm dualis transformation    ", col = "hotpink2")
    qqline(dd_log)
    
    # square root transformed
    qqnorm(dd_sq, main="Blood alcohol content of drunk drivers\nqq plot - Square root transformation", col = "hotpink2")
    qqline(dd_sq)
    
    # reciprocal transformed
    qqnorm(dd_rt, main="Blood alcohol content of drunk drivers\nqq plot - Reciprocal transformation", col = "hotpink2")
    qqline(dd_rt)
    
    par(mfrow = c(1,1)) # restore parameters
  dev.off()
  
  # Looking at the plots, which transformation looks best?
  # Looking at the plots, the histogram of the reciprocal transformed data looks best, as it is the most symmetrical of them. The others are more skewed (to the left) and have a more triangular shape, with a steep slope on the left side. The reciprocal transformed data is mre symmetrical, although it is also skewed, this time to the right, yet not as severly as the others. The slopes are also more uniform and not as steep.
  # Looking at the qq plots, the logarith dualis transformed data looks marginally better, as it fits closer to the identity line. At both the upper and lower end, the data points of all transformations veer off the identity line.
  # We conclude that none of these transformations looks ideal, but the reciprocal transformed data looks like the closest approximation to normally distributed data.

  

# 2. Two sample problem [R+]
# During a study, the weight of 458 mice (229 male, 229 female) was measured [2] and compared. The corresponding dataset "mice_study.csv" (mice_male and mice_female) can be downloaded from the course homepage. The weights of both groups were measured in grams.
# Assume a normal distribution.
  
  # import data
  mice_study <- read.csv("mice_study.csv", row.names = 1)

# [1P] a.
# Use hypothesis testing (formulas for the test statistic and p-value are given in the formula sheet) to determine if there is a difference in the mean weight between both sexes for sigma^2 x != sigma^2 y.
# State the hypotheses and compute the appropriate t-test with alpha = 0.05. Do not use a predefined R package. This task should be done on paper (it is allowed to use R as a calculator) using a t-table. State your results and do not forget to refer to the given significance level.
  
  # do this on paper!
  
  # Hypotheses
    # H0: The mean weight of the sexes is equal.
        # x_m = x_f
    # H1: The mean weight of the sexes is not equal.
        # x_m != x_f
  
  # Calculation
  alpha <- 0.05
  n_m <- length(mice_study$mice_male)
  n_f <- length(mice_study$mice_female)
  x_m <- mean(mice_study$mice_male)
  x_f <- mean(mice_study$mice_female)
  s_m <- sd(mice_study$mice_male)
  s_f <- sd(mice_study$mice_female)
  
  t <- (x_m - x_f)/sqrt((s_m^2/n_m) + (s_f^2/n_f))
  
  df <- ((((s_m^2)/n_m) + ((s_f^2)/n_f))^2) / (((1/(n_m-1))*((s_m^2)/n_m)^2) + ((1/(n_f-1))*((s_f^2)/n_f)^2))
  
    # t=8.729956
    # df=423.747848
  
  
  # Look at the table:
  # df=500, p=0.975
  # t_(500, 1-(0.05/2))=1.97
  
  # K = (-infinity : -1.97)U(1.97 : infinity)
  # t is an element of K
  
  # We reject the H0! We can say with 95% confidence that the means of the weights of the sexes are not equal.
  
  
# [1.5P] b.
# Use hypothesis testing functions in R (in this case t.test()) to check your result of 2a. and to determine if the mean weight of male mice is i) greater than or ii) less than the one of female mice for sigma^2 x != sigma^2 y. State your results based on the three types of hypotheses and draw a conclusion.
  
  # Summary
  summary(mice_study)
  
  # Is there a difference in mean weights?
  # Two-sided independent samples t-test
  
  # Hypothesis
  # H0: x_m = x_f
  # H1: x_m != x_f
  
  # t-test
  t.test(mice_study$mice_male, mice_study$mice_female, alternative = "two.sided", var.equal = FALSE)
  
          # Welch Two Sample t-test
          # 
          # data:  mice_study$mice_male and mice_study$mice_female
          # t = 8.7296, df = 423.75, p-value < 2.2e-16
          # alternative hypothesis: true difference in means is not equal to 0
          # 95 percent confidence interval:
          #   4.001196 6.326651
          # sample estimates:
          #   mean of x mean of y 
          # 46.54378  41.37985
  
  # Conclusion
  # The H0 that the means are equal is rejected (p-value < 2.2e-16). The means are different.
  
  # Is the mean weight of male mice greater than female mice?
  # One-sided independent samples t-test
  
  # Hypothesis
  # H0: x_m = x_f
  # H1: x_m > x_f
  
  # t-test
  t.test(mice_study$mice_male, mice_study$mice_female, alternative = "greater", var.equal = FALSE)
  
        # Welch Two Sample t-test
        # 
        # data:  mice_study$mice_male and mice_study$mice_female
        # t = 8.7296, df = 423.75, p-value < 2.2e-16
        # alternative hypothesis: true difference in means is greater than 0
        # 95 percent confidence interval:
        #   4.188787      Inf
        # sample estimates:
        #   mean of x mean of y
        # 46.54378  41.37985
  
  # Conclusion
  # The H0 that the means are equal is rejected (p-value < 2.2e-16). The males are heavier than the females.
  
  # Is the mean weight of male mice less than female mice?
  # One-sided independent samples t-test
  
  # Hypothesis
  # H0: x_m = x_f
  # H1: x_m < x_f
  
  # t-test
  t.test(mice_study$mice_male, mice_study$mice_female, alternative = "less", var.equal = FALSE)
  
          # Welch Two Sample t-test
          # 
          # data:  mice_study$mice_male and mice_study$mice_female
          # t = 8.7296, df = 423.75, p-value = 1
          # alternative hypothesis: true difference in means is less than 0
          # 95 percent confidence interval:
          #   -Inf 6.13906
          # sample estimates:
          #   mean of x mean of y 
          # 46.54378  41.37985 
  
  # Conclusion
  # The H0 that the means are equal is not rejected (p-value = 1). We cannot claim which group is heavier, we only failed to prove the mean weight of males to be less than the mean wight.

  # Final conclusion
  
  # From the three tests, we conclude that the mean weight of male mice is greater than the mean weight of the female mice.
  
  
  
# [0.5P] c.
# You would like to check if the male mice have gained weight after three months on a special diet. Therefore, the weights of the mice are measured on day one (data x) and the same group three months later (data y). Which hypothesis test would you choose and why? What is the advantage of this test? State suitable hypotheses and give the formula of the test statistics and the decision criteria!
  
  # Test used
  # Data x and y originate from the same subjects (same mice is measured on day one and three months later). Therefore, we would use the Paired t-Test.
  
  # Advantage of the paired t-test
  # The advantage of the paired t-test is that the variation between subjects in a population is eliminated, since the same subject is measured twice, once for each data set (on day one and three months later).
  
  # Hypotheses
  # Two-sided - Is the mean difference different from 0?
  # H0: Mean weight of male mice before and after the diet is the same
      # x_before = x_after
  # H1: Mean weight of male mice before and after the diet is not the same
      # x_before != x_after
  
  # One-sided - Is the mean difference greater than 0?
  # H0: Mean weight of male mice before and after the diet is the same
      # x_before = x_after
  # H1: Mean weight of male mice before and after the diet is not the same
      # x_before > x_after
  
  # One-sided - Is the mean difference less than 0?
  # H0: Mean weight of male mice before and after the diet is the same
      # x_before = x_after
  # H1: Mean weight of male mice before and after the diet is not the same
      # x_before < x_after
  
  # Formula of the test statistics
  # t=(x_d - mu_0)/((s_d)/(sqrt(n)))
  
  # Decision criteria
  # Reject H0 if p < alpha or if |t| > critical value associated with alpha.
  
  # We use it for continuous, normally distributed data that doesn't have major outliers. We use it when we have two independent measurements of the dependent variable, for example, weight of the same subject at two different time points, heart rate before and after exercise etc.
  
  # Test in R - two sided paired t-test
  t.test(data_x, data_y, alternative = "two.sided", paired = TRUE)
  
  
#   3. One sample problem [R]
# The tumor size of patients was measured after four months of treatment with a new drug. The dataset "tumor_size.csv" shows the size in cm and can be downloaded from the course homepage. The data is outlier adjusted. Assume a normal distribution. Further assume that the true mean of the untreated tumor-sizes is 6.2 cm.
  
  # import data
  tumor_size <- read.csv("tumor_size.csv", row.names = 1)

# [1P] a.
# Use hypothesis testing functions in R (in this case t.test()) to determine if the mean size of the tumors from patients treated with the new drug is i) different from or ii) lower than the given reference tumor size of untreated patients. State the hypotheses for i) and ii) and compute the appropriate t-test with alpha = 0.05.
  
  # mu0 = 6.2
  
  # Is the mean size of treated tumors different from the size of non-treated tumors?
  # Two-sided one sample t-test
  
  # Hypothesis
  # H0: mu_T = mu_N
  # H1: mu_T != mu_N
  
  # t-test
  t.test(tumor_size$tumor_size_treatment, alternative = "two.sided", mu = 6.2)
  
          # 	One Sample t-test
          # 
          # data:  tumor_size$tumor_size_treatment
          # t = -1.6878, df = 209, p-value = 0.09294
          # alternative hypothesis: true mean is not equal to 6.2
          # 95 percent confidence interval:
          #   5.702699 6.238541
          # sample estimates:
          #   mean of x 
          # 5.97062 
  
  # Conclusion
  # p-value = 0.09294 is greater then the alpha we set, therefore, we cannot reject the H0.
  
  # Is the mean size of treated tumors lower than the size of non-treated tumors?
  # One-sided one sample t-test
  
  # Hypothesis
  # H0: mu_T = mu_N
  # H1: mu_T < mu_N
  
  # t-test
  t.test(tumor_size$tumor_size_treatment, alternative = "less", mu = 6.2)

          # One Sample t-test
          # 
          # data:  tumor_size$tumor_size_treatment
          # t = -1.6878, df = 209, p-value = 0.04647
          # alternative hypothesis: true mean is less than 6.2
          # 95 percent confidence interval:
          #   -Inf 6.19516
          # sample estimates:
          #   mean of x 
          # 5.97062 
  
  # Conclusion
  # p-value = 0.04647 is lesser then the alpha we set, therefore, we reject the H0.
  
  # Final conclusion
  # The two-sided t-test is much more applicable than the one-sided t-test, which is seldom appropriate to use. With the two-sided t-test, we failed to disprove the H0 that the sample mean of treated tumors is equal to the mean of untreated tumors. However, with the one-sided test, where we asked if the tumor size decreased, we could with 95% significance reject the H0. We would still abide by the results of the two-sided t-test and advise that further tests on the treatment would be needed to conclude on its effectiveness.
  
  
  
# [1.5P] b.
# Determine a two sided 95%-confidence interval for the population mean mu with unknown variance for the treated patients. For calculating the confidence interval, which assumption has to be made? Which conclusion based on the confidence interval can you make regarding the true population mean?
  
  # mu0 = 6.2
  
  # Assumption
  # We assume our data is random, normally distributed and independent.
  
  # two-sided 95%-confidence interval for the population mean
    # mean of sample
    x <- mean(tumor_size$tumor_size_treatment)
    # sample variance
    s <- sd(tumor_size$tumor_size_treatment)
    # sample size
    n <- length(tumor_size$tumor_size_treatment)
    # degrees of freedom
    df <- length(tumor_size$tumor_size_treatment)-1
    # alpha
    alpha <- 0.05
    # t-value
    t <- abs(qt(p=(alpha/2), df=df, lower.tail = TRUE))
    # confidence interval
    ci_lower <- x - (t*(s/sqrt(n))); ci_lower
    ci_upper <- x + (t*(s/sqrt(n))); ci_upper
    
    # CI: [5.70; 6.24]
    
  # Conclusion
  # We conclude with 95% confidence that the true population mean of the tumor sizes of treated patients falls within the interval [5.70; 6.24].
  
  
#   [1P*] c.
# Bonus task:
#   What is the difference between the normal distribution and the t-distribution (Students t-distribution)?

    # The t-distribution is defined by the degrees of freedom. It is similar to the standard normal distribution, being symmetric and bell-shaped, but has heavier tails than the standard normal distribution. That is to adjust for the impact that outlier values have on the results in the smaller sample sizes we work with. At higher sample sizes and thus degrees of freedom, it gets more and more alike the normal distribution, with being identical at infinite degrees of freedom. Most consider the cut-off of when they're similar enough to approximate it to the normal distribution to be df=30 or df=50.


# 4. Multiple Testing [R]
# In a genome wide study SNPs (single nucleotide polymorphisms) associated with different psychiatric disorders (e.g., schizophrenia, bipolar disorder) have been identified [3]. The file "SNP_psychiatric.csv" can be downloaded on the course homepage. It consists of 14 occurring SNPs and their p-values.

  # import data
  snp_psych <- read.csv("SNP_psychiatric.csv")

# [1.5P] a.
# Correct the p-values using the Bonferroni correction as well as using the Benjamini-Hochberg correction. To do so, use the R function p.adjust() with the correct parameters. State a list of SNPs with a p-value <= 0.01, for the uncorrected p-values and the corrected p-values, respectively.
  
  # Bonferroni correction
  snp_psych$bonferroni_p <- p.adjust(snp_psych$p_value, method = "bonferroni", n=length(snp_psych$p_value))
  
  # Benjamini-Hochberg correction
  snp_psych$benj_hoch_p <- p.adjust(snp_psych$p_value, method = "BH", n=length(snp_psych$p_value))
  
  # List of SNPs with p-value <= 0.01
    # uncorrected p-value
    cat("SNPs with p-value <= 0.01 - uncorrected p-values")
    snp_psych$SNP[which(snp_psych$p_value <= 0.01)]
  
    # Bonferroni correction
    cat("SNPs with p-value <= 0.01 - Bonferroni corrected p-values")
    snp_psych$SNP[which(snp_psych$bonferroni_p <= 0.01)]
    
    # Benjamini-Hochberg correction
    cat("SNPs with p-value <= 0.01 - Benjamini-Hochberg corrected p-values")
    snp_psych$SNP[which(snp_psych$benj_hoch_p <= 0.01)]

# [0.5P] b.
# Which method of multiple testing is more stringent? Which type of error would you risk making if you did not correct the p-values?
  
  # More stringent
  # Bonferroni method
  
  # Which type of error would you risk without correction?
  # Type I errors (false positives)