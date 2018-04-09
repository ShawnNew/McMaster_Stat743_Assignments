alpha <- 0.05
# Read data from question as x
x <- c(18.5, 3.5, 5.3, 3.3, 9.4, 10.9, 21.6, 53.4, 5.3, 33.7)
len <- length(x)

# Exact Confident Interval
CI_exact <- c(qgamma(alpha/2, len) / (len * mean(x)), qgamma(1 - alpha/2, len) / (len * mean(x)))
# CI from LRT
#CI_LRT <- c()
# CI from Wald Tests
CI_Wald <- c((1 - qnorm(1 - alpha/2)/sqrt(len)) / mean(x), (1 + qnorm(1 - alpha/2)/sqrt(len)) / mean(x))
# CI from Score Tests
CI_Score <- c((1 - qnorm(1 - alpha/2)/sqrt(len)) / mean(x), (1 + qnorm(1 - alpha/2)/sqrt(len)) / mean(x))