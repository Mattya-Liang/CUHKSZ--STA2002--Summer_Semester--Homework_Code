# The R programming code for the Homework 1 in the course STA2002 - 2024 summer
# Q1(a):
times_test_vehicle <- c(1.75,1.84,2.12,1.92,2.62,2.35,3.09,3.15,
                        2.53,1.91,3.25,2.83)
summary(times_test_vehicle)
mean(times_test_vehicle)
var(times_test_vehicle)
sd(times_test_vehicle)
# Q1(b):
gemo_boxplot <- boxplot(times_test_vehicle, 
                        horizontal = TRUE, 
                        col = "lightblue", 
                        main = "Boxplot of the times to test a vehicle")

# Q2 Preface:
getwd() # 查询R语言的当前工作目录，因为smiles.txt文件在桌面上，所以需要将工作目录设置为桌面
setwd("C:/Users/12978/Desktop") # 将工作目录设置为桌面
smiles <- read.table("smiles.txt", header = TRUE) 
# 先找到smiles.txt文件的路径，然后将其读入R中
# 打开文件smiles.txt
names(smiles) = c("groups", "scores") # 为数据框添加列名
attach(smiles) # 将数据框smiles附加到搜索路径上

# Q2(Part 1):
# Construct histograms and stem-and-leaf plots for each of the four categories. 
# Comment on the shape of the distribution of the observations in each of the four categories. 
# Interpret and compare.

# 我们先进行直方图的绘制

tapply(scores,groups,stem) # stem-and-leaf plot

splitgroup = split(scores,groups) # split the data by groups
attach(splitgroup) # attach the split data to the search path
par(mfrow=c(2,2)) # set the layout of the plots
hist(false,main="") # plot the histogram of the false group
hist(felt,main="") # plot the histogram of the felt group
hist(miserable,main="") # plot the histogram of the miserable group
hist(neutral,main="") # plot the histogram of the neutral group

# 现在我们来绘制茎叶图

par(mfrow=c(2,2)) # set the layout of the plots
stem(false) # plot the stem-and-leaf plot of the false group
stem(felt) # plot the stem-and-leaf plot of the felt group
stem(miserable) # plot the stem-and-leaf plot of the miserable group
stem(neutral) # plot the stem-and-leaf plot of the neutral group

# Q2(Part 2):
# Obtain the 5-numerical summaries and the corresponding boxplots for all four categories. 
# Interpret and compare.

# 我们先进行5数总结的计算
# 只是简单的进行重复就可以了

summary(false)
var(false)
sd(false)

summary(felt)
var(felt)
sd(felt)

summary(miserable)
var(miserable)
sd(miserable)

summary(neutral)
var(neutral)
sd(neutral)

#现在我们进行图表的绘制，存在两种方法，一种是直接绘制，一种是使用boxplot函数
# 我们先进行直接绘制

par(mfrow=c(2,2))
boxplot(false)
title("false")
boxplot(felt)
title("felt")
boxplot(miserable)
title("miserable")
boxplot(neutral)
title("neutral")

# 或者我们使用boxplot函数去绘制

par(mfrow=c(1,1))
boxplot(scores~groups)

# Q2(Part 3):
# Based on the descriptive statistics and the different graphical displays, 
# Summarize your findings for the data analysis in the context of the problem.

# 我们现在根据数据去回答前面的两个问题：
#（a）微笑真的会增加宽大处理吗？（b）不同类型的微笑是否具有不同的效果？
# 从直方图，茎叶图，箱线图中我们可以看到，不同类型的微笑对于增加宽大处理的效果是不同的
# 通过均值和方差的计算，三种不同微笑的均值都比中性的表情要大，同时数据之间的方差差异并不大
# 由此我们判断，微笑的会增加宽大处理的可能性，不同的微笑类型对于增加宽大处理的效果是不同的

# Q3:
help("prop.table")
help("qnorm")

n1 <- 200 # CUHKSZ Student
X1 <- 150 # CUHKSZ Mac
n2 <- 250 # SUSTech Student
X2 <- 185 # SUSTech Mac

p1 <- X1 / n1 # The proportion of CUHKSZ students who own a Mac
p2 <- X2 / n2 # The proportion of SUSTech students who own a Mac

diff_proportion <-  p1 - p2

SE_diff <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
# We use SE_diff to represent the standard error of the difference between two sample proportions
# Also be the standard deviation of the sampling distribution of the difference between two sample proportions

print(paste("CUHKSZ Ratio p1 =", p1))
print(paste("SUSTech Ratio p2 =", p2))
print(paste("Sample proportion distribution =", diff_proportion))
print(paste("Sample proportion distribution Error SE =", SE_diff))
print(paste("The Normal distribution of the difference between two sample proportions is Mean",
            diff_proportion, "and Variance", SE_diff^2))
plot(seq(-0.2, 0.2, 0.001), dnorm(seq(-0.2, 0.2, 0.001), mean = diff_proportion, sd = SE_diff), 
     type = "l", col = "blue", xlab = "Difference in Proportions", ylab = "Density", 
     main = "Normal Distribution of the Difference in Proportions", lwd = 2)

# Q4(a):
# Easy to find that the sample mean and population are not always the same
# 容易发现，我们在unbiased的估计中，我们使用的是样本均值，而不是总体均值
# 并且两个变量相互独立，因此在计算estimator的时候，我们可以直接使用sample mean去代替population mean
# 也就是 miu_1 - miu_2 = Xbar_1 - Xbar_2
# 现在我们来进行计算Standard Deviation

mu1 <- Xbar1
mu2 <- Xbar2
mu_hat_diff <- mu1_hat - mu2_hat

standard_deviation <- sqrt((sigma1 ^ 2) / n1 + (sigma2 ^ 2 / n2)) 

# Q4(b):
# Now we have the biased sample and we start to calculate
# The bias of the estimator X̄1 - X̄2 for μ1 - μ2 can be calculated by taking the expected value of the estimator and subtracting the true parameter value:
# Bias(X̄1 - X̄2) = E(X̄1 - X̄2) - (μ1 - μ2) = μ1 - μ2 - (μ1 - μ2) = 0.
# The bias of this estimator is zero, indicating that it is an unbiased estimator for the difference of population means. 
# As the sample sizes n1 and n2 increase to infinity, the bias remains zero.

# Q4(c):
# To show that Sp^2 = ((n1 - 1)S1^2 + (n2 - 1)S2^2) / (n1 + n2 - 2) is an unbiased estimator of σ^2, we need to demonstrate that its expected value is equal to σ^2:
# E(Sp^2) = E(((n1 - 1)S1^2 + (n2 - 1)S2^2) / (n1 + n2 - 2)).
# Expanding the expression:
# E(Sp^2) = (n1 - 1)E(S1^2) + (n2 - 1)E(S2^2) / (n1 + n2 - 2).
# Since S1^2 and S2^2 are unbiased estimators of σ1^2 and σ2^2, respectively, their expected values are equal to the corresponding population variances:
# E(S1^2) = σ1^2 and E(S2^2) = σ2^2.
# Substituting these values:
# E(Sp^2) = (n1 - 1)σ1^2 + (n2 - 1)σ2^2 / (n1 + n2 - 2).
# Since we assume that both populations have the same variance (σ1^2 = σ2^2 = σ^2), we can simplify further:
# E(Sp^2) = (n1 - 1)σ^2 + (n2 - 1)σ^2 / (n1 + n2 - 2).
# Combining terms:
# E(Sp^2) = ((n1 - 1)σ^2 + (n2 - 1)σ^2) / (n1 + n2 - 2).
# Factoring out σ^2:
# E(Sp^2) = (n1 - 1 + n2 - 1)σ^2 / (n1 + n2 - 2).
# Simplifying:
# E(Sp^2) = (n1 + n2 - 2)σ^2 / (n1 + n2 - 2).
#Canceling out the common term (n1 + n2 - 2):
# E(Sp^2) = σ^2.
# Therefore, Sp^2 is an unbiased estimator of σ^2.

# Q5: Maximum Likelihood Estimation

# Q5(a): Using MOM to estimate the parameter of the exponential distribution
# 1. MOM: Mean of the exponential distribution is 1 / lambda
# 2. We find that the lambda hat is 1 / Xbar

# Q5(b): Using MLE to estimate the parameter of the exponential distribution
# 1. MLE: The likelihood function of the exponential distribution is lambda * exp(-lambda * x)
# 2. We use the formula to calculate the MLE of lambda
# 3. The MLE of lambda is the all multipicaltion of the function f(x) = lambda * exp(-lambda * x)
# 4. The l(x) = log[L(x)], then we use the partial derrivative to calculate the MLE of lambda
# 5. Let the partial derrivative of l(x) to lambda be 0, then we can get the MLE of lambda
# 6. The MLE of lambda is 1 / Xbar

# Q5(c):
Sample <- c(3.8,3.24,1.4,1.22,4.5,4.6)
mean_sample <- mean(Sample)
Xbar <- 1 / mean_sample
cat("The MOM of lambda is", Xbar, "\n")
cat("The MLE of lambda is", Xbar, "\n")

# Q6:(a)

# Q6:(b)

# Q7:
