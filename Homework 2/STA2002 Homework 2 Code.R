# The R programming code for the Homework 2 in the course STA2002 - 2024 summer

# Q1(a):Calculate the bias of each point estimate. Is any one of them unbiased?

# Easy to calculate the bias of each point estimate by subtracting the true value from the point estimate.
# The bias of each point estimate is calculated as follows:
# miu_hat1 = 1/3*(X1+X2+X3) = 1/3*(miu+miu+miu) = miu
# miu_hat2 = 1/4*X1 + 1/3*X2 + 1/5*X3 = 47/60*miu
# Bias of miu_hat1 = miu - miu = 0
# Bias of miu_hat2 = 47/60*miu - miu = -13/60*miu
# Therefore, miu_hat1 is unbiased, while miu_hat2 is biased.

# Q1(b):Calculate the variance of each point estimate. Which point estimate has the smallest variance?
# Var(X1) = 7, Var(X2) = 13, and Var(X3) = 20.
# The variance of each point estimate is calculated as follows:
# Var(miu_hat1) = 1/3^2*Var(X1) = 40/9
# Var(miu_hat2) = 1/4^2*Var(X1) + 1/3^2*Var(X2) + 1/5^2*Var(X3)
# = 1/16*7 + 1/9*13 + 1/25*20 = 1.6
# Therefore, miu_hat1 has the smallest variance.

# Q1(c):Calculate the mean squared error of each point estimate. 
# Which point estimate has the smallest mean squared error?
# The mean squared error of each point estimate is calculated as follows:
# MSE(miu_hat1) = Var(miu_hat1) + Bias(miu_hat1)^2 = 40/9
# MSE(miu_hat2) = Var(miu_hat2) + Bias(miu_hat2)^2 = 1.6 + (-13/60*miu)^2
# Therefore, miu_hat1 has the smallest mean squared error.

# Q2:Find the maximum likelihood estimator of the unknown parameter θ
# We have X1,X2,...,Xn ~ Exp(θ)
# Where the density function of the distribution is f(x|θ) = exp(-(x-θ)) for x >= 0 and 0 otherwise
# The likelihood function is L(θ) = exp(-n*θ)*exp(-sum(Xi))
# The log-likelihood function is l(θ) = -n*θ - sum(Xi)
# The maximum likelihood estimator of the unknown parameter θ is the value of θ that maximizes the log-likelihood function
# To find the maximum likelihood estimator of θ, we take the derivative of the log-likelihood function with respect to θ and set it equal to zero
# The derivative of the log-likelihood function with respect to θ is d/dθ(l(θ)) = -n + sum(Xi)
# Setting the derivative equal to zero, we get -n + sum(Xi) = 0
# Which only exists when n = 0, which is impossible
# Thus, we consider again the function again
# The above function is maximum when exp(θn) is maximum, but θ is a positive number
# Hence we can say that the θ_hat = min(X1,X2,...,Xn) is the maximum likelihood estimator of the unknown parameter θ


# Q3:Gaussian Mixture Model(GMM)
# The Gaussian Mixture Model (GMM) is a probabilistic model that assumes that the data is generated from a mixture of several Gaussian distributions.
# The GMM is defined by the following parameters:
# 1. The number of components (K) in the mixture model.
# 2. The mean (μ) and variance (σ^2) of each component.
# 3. The mixing coefficients (π) that determine the probability of each component being selected.

# Q3(a):Derive the joint density of (X;K). State clearly the support of (X;K) in the joint density. 
# Hint: consider conditional distribution and the law of total probability.
# The joint density of (X;K) can be derived as follows:
# 我们知道：f(X,K)可以被写成f(X = x,K = k)的形式
# 而这个形式又可以通过条件概率的形式写成f(X = x|K = k)*f(K = k)的形式
# f(X = x|K = k)是给定K的条件下X的密度，这是一个混合了K个高斯分布的密度
# f(K = k)是K的密度，这是一个离散的分布，取值为0或1
# 于是我们可以使用一个条件函数也就是 f_{X|K}(x|k) pi_k来表示
# 此处的pi_k表示的是K是k的概率，也就是所谓的pi_0或者pi_1
# 因此，我们可以表示联合函数为f(X,K) = f_{X|K}(x|k) pi_k
# 也就是N(x|mu_k, sigma_k^2)*pi_k


# Q3(b):
# The likelihood function for the i.i.d. sample of size n is given by:
# L(pi_0, mu_0, sigma_0^2, mu_1, sigma_1^2)

# Q4:
# Data: 15.2 14.2 14.3 14.2 14.0 13.5 12.2 11.8 14.4 12.5 15.2
# Assume the data are random samples from a normal distribution, with the standard deviation known to be = 0.5.
# In this case, the sample mean is an unbiased estimator of the population mean, and the sample variance is an unbiased estimator of the population variance.
# In other words, we only know the variance of the population, and unknown population mean
# The sample mean is calculated as follows:
# X_bar = (15.2 + 14.2 + 14.3 + 14.2 + 14.0 + 13.5 + 12.2 + 11.8 + 14.4 + 12.5 + 15.2)/11 = 13.77
# The sample variance is calculated as follows:
# S^2 = ((15.2-13.7)^2 + (14.2-13.7)^2 + (14.3-13.7)^2 + (14.2-13.7)^2 + (14.0-13.7)^2 + (13.5-13.7)^2 + (12.2-13.7)^2 + (11.8-13.7)^2 + (14.4-13.7)^2 + (12.5-13.7)^2 + (15.2-13.7)^2)/(11-1) = 0.5
# The standard deviation is known to be 0.5.
# Therefore, the sample mean is 13.7, the sample variance is 0.5, and the standard deviation is 0.5.

# Q4(a): Construct a 99% two-sided confidence interval on the mean temperature.
# The confidence interval is calculated as follows:
# CI = X_bar +/- Z*sigma/sqrt(n)
# We can easily use: X_bar = 13.77, alpha = 0.01, sigma = 0.5 to calculate
# Thus, We can get the CI: 14.87 >= miu >= 12.67

# Q4(b): Construct a 95% lower-confidence bound on the mean temperature.
# The confidence interval is calculated as follows:
# CI = X_bar + 
# Thus, We can get the lower bounded CI: positive infinite >= miu >= 12.49

# Q4(c):Suppose that we wanted to be 95% confident that the error in estimating the mean temperature is less than 2 degrees Celsius. What sample size should be used?
# We just want the minimun size of the sample
# The length of CI is given by 2*Z*sigma/sqrt(n)

# Now we try to use R to calculate the answer:
# We will use t distribution function to test, we first install the package of "pwr"
install.packages("pwr")
library(pwr)

# For Question (a):
Temperatures <- c(15.2,14.2,14.3,14.2,14.0,13.5,12.2,11.8,14.4,12.5,15.2) 
# Input the data of temperature
CI_a <- t.test(Temperatures, conf.level = 0.99)$conf.int
CI_a

# For Question (b):
lower_CI <- t.test(Temperatures, conf.level = 0.95, alternative = "less")$conf.int
lower_CI

# For Question (c):
required_n <- pwr.t.test(d = 2/0.5, sig.level = 0.05, power = 0.95, type = "two.sample")$n
required_n

# Q5:
# I will just show the code 

ct_scans <- c(2.31, 2.09, 2.36, 1.95, 1.98, 2.25, 2.16, 2.07, 1.88, 1.94, 1.97, 2.02)
s <- sd(ct_scans)
s

n <- length(ct_scans) # Sample size

df <- n - 1 # Degrees of freedom

conf_level <- 0.95 # Confidence level with two-sided at significance of 95%

chi_lower <- qchisq(1 - conf_level/2, df)
chi_upper <- qchisq(conf_level/2, df)

# Critical values of chi-square

CI_sd <- c(sqrt((df * s^2) / chi_upper), sqrt((df * s^2) / chi_lower))
CI_sd

  
