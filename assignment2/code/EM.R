# R code to implement the EM algorithm for the two examples 
# done in class on February 1, 2018
#
ll <- function(p, x) {
  sum(log(p*dnorm(x, mean=10, sd=1) + (1-p)*dnorm(x, mean=13, sd=1)))
}

EM <- function(x, init, ll, eps=1e-12, maxit=100) {
  # A function to implement the EM algorithm for the censored
  # exponential example. The function takes the observed data, 
  # an intital estimate, and the log likelihood function as well
  # as optional tolerance and maximum iteration parameters.
  out <- matrix(NA, nrow=maxit+1, ncol=4)
  out[1,1:2] <- c(init, ll(init, x))
  i <- 1
  continue <- TRUE
  old <- init
  n <- length(x)
  while (continue) {
    new <- sum(old*dnorm(x, mean=10, sd=1) / (old*dnorm(x, mean=10, sd=1) + (1-old)*dnorm(x, mean=13, sd=1))) / 20
    out[i+1,1:2] <- c(new, ll(new,x))
    out[i+1,3] <- abs(new-old)
    out[i+1,4] <- out[i+1,2]-out[i,2]
    i <- i+1
    old <- new
    continue <- (i<=maxit & out[i,4]>eps)
  }
  out <- out[1:i,]
  return(list(est=new, trace=out))
}


x <- c(9.29, 12.86, 9.73, 11.45, 10.13, 9.55, 9.00, 9.78,
       12.74, 9.49, 9.70, 13.38, 9.08, 13.35, 9.33, 14.31,
       10.10, 10.03, 10.86, 9.76) #  the given dataset
p.mle <- EM(x, 0.5, ll)
