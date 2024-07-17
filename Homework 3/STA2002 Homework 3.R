# The R programming code for the Homework 3 in the course STA2002 - 2024 summer
# Q1：
# Recent Gallup Poll estimates that 88% Americans believe that cloning humans is morally unacceptable. 
# Results are based on telephone interviews with a randomly selected national sample of n = 1000 adults, aged 18 and older, conducted May 2-4, 2004.

# Q1(a):
# Find 95% condence interval for the true proportion? Does 0.9 fall in the interval?
# The 95% confidence interval for the true proportion is given by:
# p̂ ± z*sqrt(p̂(1-p̂)/n)
# where p̂ is the sample proportion, z is the z-score corresponding to the desired confidence level, and n is the sample size.
# The z-score for a 95% confidence level is 1.96.
# The sample proportion is 0.88 and the sample size is 1000.
# Therefore, the 95% confidence interval for the true proportion is given by:
# 0.88 ± 1.96*sqrt(0.88*(1-0.88)/1000)
# 0.88 ± 1.96*sqrt(0.88*0.12/1000)  
# 0.88 ± 1.96*sqrt(0.1056/1000)
# 0.88 ± 1.96*sqrt(0.0001056)
# 0.88 ± 1.96*0.010273
# 0.88 ± 0.0201
# 0.8599 to 0.9001
# Therefore, the 95% confidence interval for the true proportion is 0.8599 to 0.9001.
# Since 0.9 falls within this interval, we can conclude that the true proportion of Americans who believe that cloning humans is morally unacceptable is likely to be between 0.8599 and 0.9001.
# The answer is yes.

# 给定的数据和参数
p_hat <- 0.88  # 样本比例
n <- 1000      # 样本大小
alpha <- 0.05  # 显著性水平 (1 - 置信水平)

# 计算标准误差
SE <- sqrt(p_hat * (1 - p_hat) / n)

# 计算临界值（Z分位数）
z_critical <- qnorm(1 - alpha/2)

# 计算置信区间的边界
margin_of_error <- z_critical * SE

# 计算置信区间
lower_bound <- p_hat - margin_of_error
upper_bound <- p_hat + margin_of_error

# 输出结果
cat("95%置信区间为: (", lower_bound, ",", upper_bound, ")\n")

# 检查0.9是否在置信区间内
if (lower_bound <= 0.9 && upper_bound >= 0.9) {
  cat("0.9 在置信区间内.\n")
} else {
  cat("0.9 不在置信区间内.\n")
}
# Q1(b):Pretend that you want to replicate Gallup's inquiry in the Shenzhen area. 
# What sample size is needed so that the length of a 95% confidence interval for the unknown proportion of people in the area who believe that cloning humans is morally unacceptable does not exceed 0.02?
# The formula for the margin of error in a confidence interval is given by:
# z*sqrt(p̂(1-p̂)/n)
# where z is the z-score corresponding to the desired confidence level, p̂ is the sample proportion, and n is the sample size.
# In this case, we want the margin of error to be less than or equal to 0.02. Therefore, we can set up the following inequality:
# z*sqrt(p̂(1-p̂)/n) ≤ 0.02
# To solve for n, we can rearrange the inequality as follows:
# n ≥ z^2 * p̂(1-p̂) / (0.02)^2
# Given that z = 1.96 for a 95% confidence level, p̂ = 0.88, and we want the margin of error to be 0.02, we can plug in these values to calculate the required sample size.
# n ≥ 1.96^2 * 0.88(1-0.88) / 0.02^2
# n ≥ 3.8416 * 0.88(0.12) / 0.0004
# n ≥ 3.381408 / 0.0004
# n ≥ 8453.52
# Therefore, the sample size needed so that the length of a 95% confidence interval for the unknown proportion of people in the Shenzhen area who believe that cloning humans is morally unacceptable does not exceed 0.02 is 8454.
# 给定参数
alpha <- 0.05  # 显著性水平 (1 - 置信水平)
margin_of_error <- 0.02  # 置信区间的最大长度

# 标准正态分布的临界值 (95% 置信水平的 Z 分位数)
z_critical <- qnorm(1 - alpha/2)

# 初始估计的比例 p （使用0.5作为保守估计）
p <- 0.88

# 计算所需的样本大小
n <- ceiling((z_critical / margin_of_error)^2 * p * (1 - p))

# 输出结果
cat("所需的样本大小为:", n, "\n")

# Q2:
# A manufacturer is interested in the output voltage of a power supply used in a PC. 
# Output voltage is assumed to be normally distributed, with standard deviation 0.25 volt, and the manufacturer wishes to test H0 : miu = 5 volts against H1 : miu != 5 volts, using n = 8 units.
# Suppose we set the acceptance region to be 4.85 =< x =< 5.15.

# Q2(a): Find the Type I error (alpha) of the test.

# 我们先画一个H0的正态分布函数，并且标注出接受域

mean_value <- 5
sd_value <- 0.25

x <- seq(3, 7, length.out = 100)  # 选择适当的范围
y <- dnorm(x, mean = mean_value, sd = sd_value)

plot(x, y, type = "l", lwd = 2, col = "blue", xlab = "H0: Null Hypothesis", ylab = "Density", 
     main = "Normal Distribution with Mean 5 and SD 0.25")

legend("topright", legend = "N(5, 0.25)", col = "blue", lty = 1, lwd = 2)

segments(4.85, 0, 4.85, dnorm(4.85, mean = mean_value, sd = sd_value), col = "red", lwd = 2)
segments(5.15, 0, 5.15, dnorm(5.15, mean = mean_value, sd = sd_value), col = "red", lwd = 2)
segments(4.85, dnorm(4.85, mean = mean_value, sd = sd_value), 5.15, dnorm(5.15, mean = mean_value, sd = sd_value), col = "red", lwd = 2)

# 现在我们计算alpha的值，也就是第一类错误的概率，也就是落在拒绝域内的概率即红线以外的概率
alpha <- pnorm(4.85, mean = mean_value, sd = sd_value) + 1 - pnorm(5.15, mean = mean_value, sd = sd_value)

# 正常计算的话：
# 也就是计算x小于4.85和x大于5.15的概率之和
# 我们分成两步计算
p1 <- pnorm(4.85, mean = mean_value, sd = sd_value)
p2 <- 1 - pnorm(5.15, mean = mean_value, sd = sd_value)
alpha <- p1 + p2
cat("Type I error (alpha) of the test is:", alpha, "\n")

# Q2(b): (b) Find the Type-II error (beta) of the test for detecting a true mean output voltage of 5.1 volts.
# 这是第二类错误，即我们接受了H0，但是实际上H1是正确的，也就是我们没有拒绝H0，但是实际上H1是正确的
# 我们需要计算在H1为真的情况下，我们没有拒绝H0的概率
# 我们先画一个H1的正态分布函数，均值是5.1，而且标注出接受域
mean_value_h1 <- 5.1
sd_value_h1 <- 0.25
x <- seq(3, 7, length.out = 100)  # 选择适当的范围
y <- dnorm(x, mean = mean_value_h1, sd = sd_value_h1)
plot(x, y, type = "l", lwd = 2, col = "blue", xlab = "H1: Alternative Hypothesis", ylab = "Density", 
     main = "Normal Distribution with Mean 5.1 and SD 0.25")
legend("topright", legend = "N(5.1, 0.25)", col = "blue", lty = 1, lwd = 2)
abline(v = 4.85, col = "red", lwd = 2)
abline(v = 5.15, col = "red", lwd = 2)

# 我们要计算的第二类错误的概率，也就是红色区域的面积
# 直接计算就可以了
beta <- pnorm(5.15, mean = mean_value_h1, sd = sd_value_h1) - pnorm(4.85, mean = mean_value_h1, sd = sd_value_h1)
cat("Type-II error (beta) of the test for detecting a true mean output voltage of 5.1 volts is:", beta, "\n")

# Q2(c): Find the p-value when the observed statistic is (i) x_bar = 5.2 and (ii) x_bar = 4.7, respectively.
# 我们要计算的是p值，也就是在H0为真的情况下，观察到的统计量的概率
mean_value <- 5
sd_value <- 0.25

x <- seq(3, 7, length.out = 100)  # 选择适当的范围
y <- dnorm(x, mean = mean_value, sd = sd_value)

plot(x, y, type = "l", lwd = 2, col = "blue", xlab = "H0: Null Hypothesis", ylab = "Density", 
     main = "Normal Distribution with Mean 5 and SD 0.25")

legend("topright", legend = "N(5, 0.25)", col = "blue", lty = 1, lwd = 2)

segments(4.85, 0, 4.85, dnorm(4.85, mean = mean_value, sd = sd_value), col = "red", lwd = 2)
segments(5.15, 0, 5.15, dnorm(5.15, mean = mean_value, sd = sd_value), col = "red", lwd = 2)
segments(4.85, dnorm(4.85, mean = mean_value, sd = sd_value), 5.15, dnorm(5.15, mean = mean_value, sd = sd_value), col = "red", lwd = 2)

# 先从x_bar = 5.2开始
p_value_1 <- 2 * pnorm(5.2, mean = mean_value, sd = sd_value, lower.tail = FALSE)
cat("The p-value when the observed statistic is x_bar = 5.2 is:", p_value_1, "\n")

# 然后计算x_bar = 4.7的情况
p_value_2 <- 2 * pnorm(4.7, mean = mean_value, sd = sd_value)
cat("The p-value when the observed statistic is x_bar = 4.7 is:", p_value_2, "\n")



