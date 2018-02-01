eff <- function(n) {
  cn <- sqrt((n-1)/2) * gamma((n-1) / 2) / gamma(n/2)
  res = 1 / (2*n*(cn^2-1))
  res
}

ratio <- function(n) {
  cn <- sqrt((n-1)/2) * gamma((n-1) / 2) / gamma(n/2)
  res = (cn^2-1) / (2*(1-1/cn))
  res
}

n <- 100
effResult <- matrix(ncol=1,nrow=n)
ratioResult <- matrix(ncol=1, nrow=n)
for (i in 1:n) {
  effResult[i] <- eff(i)
  ratioResult[i] <- ratio(i)
}

plot(effResult, type = "p", xlab="n", ylab="efficiency", pch=20)
title("efficiency with respect to n")
plot(ratioResult, xlab="n", ylab="ration", pch=20)
title("ration with respect to n")
