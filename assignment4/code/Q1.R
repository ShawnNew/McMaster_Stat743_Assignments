alpha <- 0.05
# Read data from question as x
x <- c(18.5, 3.5, 5.3, 3.3, 9.4, 10.9, 21.6, 53.4, 5.3, 33.7)
len <- length(x)

# Exact Confident Interval
CI_exact <- c(qgamma(alpha/2, len) / (len * mean(x)), qgamma(1 - alpha/2, len) / (len * mean(x)))
# CI from LRT
f <- function(lambda) lambda*mean(x)-log(lambda)-(0.5*qchisq(1-alpha, 1)+len*log(mean(x))+len)/len
CI_LRT <- c(uniroot(f, c(0,0.05), tol=0.0001)$root, uniroot(f, c(0.05,0.15), tol=0.0001)$root)
# CI from Wald Tests
CI_Wald <- c((1 - qnorm(1 - alpha/2)/sqrt(len)) / mean(x), (1 + qnorm(1 - alpha/2)/sqrt(len)) / mean(x))
# CI from Score Tests
CI_Score <- c((1 - qnorm(1 - alpha/2)/sqrt(len)) / mean(x), (1 + qnorm(1 - alpha/2)/sqrt(len)) / mean(x))


## Draw line to specify the interval of root
dd <- seq(0, 0.25, by=0.000001)
tt <- f(dd)
df <- data.frame(dd,tt)
g <- ggplot(df,aes(dd,tt))
g <- g+geom_line(col='red')   #red line
g