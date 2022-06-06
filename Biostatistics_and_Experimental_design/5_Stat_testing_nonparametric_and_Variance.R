#Excercise 5

# 1. One sample problem [R+]

# During a study the blood cholesterol (in mg/dL) of men in their thirties who consume high amounts of processed foods was measured.
# cholesterol = 163; 145; 175; 147; 136; 151; 166; 162; 187
# Assume the sample to satisfy a continuous distribution, which is symmetrical around the median. By using the Wilcoxon signed-rank test for a median and a significance level of α = 0.05 determine the following points:

# [1.5P] a.
# Assess if the median cholesterol level differs significantly from m0 = 145 (the cholesterol level of a healthy male). Use the given formulas (test statistics W+ and critical region K) and the Wilcoxon-table from the homepage to solve this task. Do not use a predefined R-function. The calculations (using the given formulas), however, can be done in R or on paper. State your conclusions by referring to your calculated test statistics and the critical region.

# b) [1P] b.
# Assess if the median cholesterol level is i) less than m0 = 145 and ii) greater than m0 = 145. Solve this task by using the predefined R-function wilcox.test(). State your hypotheses and results for both cases. For your conclusions, refer to your calculated test statistics and the p-values (which are computed by the R-function).

# i) less than m0 = 145

  # Hypotheses
  # H0: m=m0 (The median cholesterol level is equal to the cholesterol level of a healthy male) 
  # H1: m!=m0 (The median cholesterol level is not equal to the cholesterol level of a healthy male) 

chol_m<- c(163, 145, 175, 147, 136, 151, 166, 162, 187)

wilcox.test(chol_m,  var.equal = FALSE)

wilcox.test(chol_m, alternative = "less", var.equal = FALSE)

# ii) greater than m0 = 145.
wilcox.test(chol_m, alternative = "greater", var.equal = FALSE)

# [0.5P] c.
# Summarize your results based on the three hypotheses tested in a. and b.
# Which conclusion can you draw regarding the measured cholesterol level of men who consume a lot of processed food in comparison to the cholesterol level of a healthy male?

#2) Two sample problem [R+]

# To assess the contribution of age on the blood cholesterol levels, the blood cholesterol of two different age groups (A = thirties, B = sixties) with similar diets was measured. Assume that the samples of the two groups are independent and satisfy a continuous distribution, which is symmetrical around the median.
# A = 160; 171; 153; 151; 175; 159; 158; 154; 165
# B = 162; 176; 183; 185; 157; 167; 164; 175; 153
# By using a non-parametric test and a significance level of α = 0.05 determine the following points:

# [1.5P] a.
# Assess if there is a difference in the median of the two age groups (A and B), using the Wilcoxon rank sum test. Solve this task by using the given formulas (critical region K and test statistic WN) and tables from the homepage and not by using a predefined R-function. The calculations can be done in R or on paper. In addition, insert a table containing the combined and sorted ranks.

# [1P] b.
# Determine if the median cholesterol level of age group A (= in their thirties) is i) greater or ii) less than the cholesterol level of age group B (= in their sixties) using the Wilcoxon rank sum test. Solve this task by using the given formulas and tables on the homepage and not by using a predefined R function. Calculations can be done in R or on paper.

# [1P] c.
# Check your results of 2a. by using the R-function wilcox.test(), which computes the Mann-Whitney-U-test statistic (not WN but you can convert using the given formulas!).

#Man Withney U - Test

  # Hypotheses
  # H0: ma=mb (The median of group A and B are equal)
  # H1: ma!=mb (The median of group A and B are not equal)

thirties<-c(160,171,153,151,175,159,158,154,165)
sixties<-c(162,176,183,185,157,167,164,175,153)

alpha<-0.05

cholesterol_lvl <- data.frame(thirties,sixties)
colnames<-c("thirties","sixties")

wilcox.test(cholesterol_lvl$thirties, cholesterol_lvl$sixties, paired=FALSE)

# [0.5P] d.
# Summarize your results from a. and b. Which conclusion can you draw regarding the median cholesterol level of the two different age groups?

#   [0.5P*] e.
# Bonus task:
#   If you could assume a normal distribution, which parametric test would you chose for this exercise and why?


# 3. Analysis of Variance [R]

# During a study, participants tried losing weight with one of three diets. Their weight at the beginning, and after 6 weeks of sticking to the diet were measured [1], and the amount of lost weight was calculated.
# D1 = 3.8; 4.0; 0.7; 2.9; 2.8; 2.0; 2.0; 8.5; 1.9; 3.1; 1.5; 3.0; 3.6; 0.9
# D2 = 2.1; 2.0; 1.7; 4.3; 6.0; 0.6; 2.7; 3.6; 3.0; 2.0; 2.2; 1.7; 3.3; 0.5
# D3 = 7.0; 5.6; 6.4; 6.8; 7.8; 8.4; 6.8; 7.2; 7.0; 7.3; 9.5; 7.6; 6.1; 6.3
# The success of the three diets should be analysed, using ANOVA. For the following tasks, assume the samples to be normally distributed.

# [2P] a.
# Investigate if the three diets show a difference in the overall amount of lost weight based on ANOVA testing. State your null and alternative hypothesis. Use the given formula to calculate the F-value. Verify your results by using the R-function aov() (for more information on the function use ?aov). Include the summary table of aov in your report (Hint: Use the command summary()). What can you conclude from the test results?

  # Hypotheses
  # H0: µ1=µ2=µ3 (The three samples come from populations with equal mean)
  # H1: µ1=µ2!=µ3; µ1!=µ2=µ3; µ1!=µ2!=µ3 (at least one mean of the three samples differs from the others) 


#Loading data
D1 <- c(3.8, 4, 0.7, 2.9, 2.8, 2, 2, 8.5, 1.9, 3.1, 1.5, 3, 3.6, 0.9)
D2 <- c(2.1, 2, 1.7, 4.3, 6, 0.6, 2.7, 3.6, 3, 2, 2.2, 1.7, 3.3, 0.5)
D3 <- c(7, 5.6, 6.4, 6.8, 7.8, 8.4, 6.8, 7.2, 7, 7.3, 8.5, 7.6, 6.1, 6.3)

#creating dataframe
data <- stack(data.frame(D1, D2, D3))
data1 <- data.frame(D1, D2, D3)
data2 <- data.frame(c(3.8, 4, 0.7, 2.9, 2.8, 2, 2, 8.5, 1.9, 3.1, 1.5, 3, 3.6, 0.9),
                    c(2.1, 2, 1.7, 4.3, 6, 0.6, 2.7, 3.6, 3, 2, 2.2, 1.7, 3.3, 0.5),
                    c(7, 5.6, 6.4, 6.8, 7.8, 8.4, 6.8, 7.2, 7, 7.3, 8.5, 7.6, 6.1, 6.3))
colnames(data2) <- c("D1", "D2", "D3")


k <- 3
n <- length(D1) +length(D2) + length(D3)

#Calculating variance within groups

var_within <- 0

for (i in 1:length(D1))
{ 
  var_within <- var_within + ((D1[i] - mean(D1))^2 + (D2[i] - mean(D2))^2 + (D3[i] - mean(D3))^2)
}

#Normalizing
var_within <- var_within / dim(data)[1]
var_within

#Calculating variances between groups
var_between <- ((length(D1) * (mean(D1) - mean(data$values))^2) + 
                 (length(D2) * (mean(D2) - mean(data$values))^2) + 
                 (length(D3) * (mean(D3) - mean(data$values))^2)) / dim(data)[1]
var_between

#Calculating F value

F_value1 <- ((n - k) * var_between) / ((k - 1) * var_within); F_value1

#Calculating ANOVA
anova <- aov(values ~ ind, data = data)
summary(anova)

# [1P] b.
# Create a boxplot diagram of all three samples (all three boxplots should be depicted in one plot). Use the identical parameter definition as with aov to infer the boxplot. What can be seen in the diagram (e.g., which sample has the highest median? What can be said with regards to the variances of the samples?)? Give the relations of the variances of the three samples, as well as the relations of the medians of the three samples in mathematical notation (e.g., x < y < z).

png("diet_boxplot.png")
boxplot(data$values ~ data$ind, col=c("lightblue", "pink", "lightgreen"), names = c("Diet 1", "Diet 2", "Diet 3"), xlab = "Diet plan", ylab = "Average weight loss [units not stated]", main = "Average weight loss by diet type")
dev.off()

# [1P*] c.
# Bonus task:
#   Use the F-test for two samples to check at a significance level α = 0.05 if the variances of diets D1 and D3 differ (two-sided F-test). Use the given formulas and the F-statistic table from the homepage to find out the critical range. Check your results in R using the function var.test() (Hint: var.test() does not give the critical range). Sum up your results in words and give the relation between the two variances in mathematical notation.

# Hypotheses
# H0: Tau� = 1, var_D1 = var_D3 (The variances of groups D1 and D3 are equal)
# H1: Tau� != 1, var_D1 != var_D3 (The variances of groups D1 and D3 are not equal)

# F-test by hand - two sided F-test
alpha=0.05

F_value <- var(D1)/var(D3)
F_value

df_D1 <- length(D1) - 1; df_D1
df_D3 <- length(D3) - 1; df_D3

F_10_10_0.975 <- 3.717

F_10_10_0.025 <- 1/F_10_10_0.975

cat("K = [0;", F_10_10_0.025, ") U (", F_10_10_0.975, "; inf)", sep = "")

var.test(D1, D3)

