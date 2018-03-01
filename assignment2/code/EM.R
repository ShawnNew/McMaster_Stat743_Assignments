# R code to implement the EM algorithm for the two examples 
#
ll <- function(p, x) {
  sum(log(p*dnorm(x, mean=10, sd=1) + (1-p)*dnorm(x, mean=13, sd=1)))
}

EM <- function(x, init, ll, eps=1e-5, maxit=100000) {
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
