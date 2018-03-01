# R code to implement the Newton-Raphson algorithm
#


U <- function(p, x)
  # The data x are ten-dimensional
  # The derivative of the log-likelihood for pi 
  sum((dnorm(x, mean=10, sd=1) - dnorm(x, mean=13, sd=1)) / 
        (p * dnorm(x, mean=10, sd=1) + (1-p) * dnorm(x, mean=13, sd=1)))

J <- function(p, x) 
  # The data x are ten-dimensional
  # The NEGATIVE of the second derivative of the log-likelihood
  # for pi
  sum((dnorm(x, mean=10, sd=1) - dnorm(x, mean=13, sd=1))**2 / 
         (p * dnorm(x, mean=10, sd=1) + (1-p) * dnorm(x, mean=13, sd=1))**2)
  #trigamma(al)-1/al

newton <- function (th0, x, U, J, eps=1e-5, maxit=100000) {
  # A general function to implement a one-dimensional 
  # Newton-Raphson algorithm to solve the likelihood equation.
  out <- matrix(NA, nrow=maxit+1, ncol=4) # Output matrix
  out[1,1:3] <- c(th0, U(th0, x), J(th0, x))
  continue <- TRUE
  iter <- 1
  while (continue) {  # While loop to iterate the algorithm
    # Get the updated estimate
    theta.new <- out[iter,1]+out[iter,2]/out[iter,3]
    iter <- iter+1
    out[iter,1:3] <- c(theta.new, U(theta.new, x), 
                       J(theta.new,x))
    out[iter, 4] <- abs(out[iter,1]-out[iter-1, 1])
    # Now check to see if convergence has been achieved.
    # We terminate if BOTH U(theta.new)<eps and 
    # |theta.new-theta.old|<eps or the max.iter iterations 
    # have been completed.
    continue <- (iter<=maxit & out[iter,4]>eps)
  }
  out <- out[1:iter,]
  return(list(est=out[iter,1], trace=out))
}

# compute and draw picture
x <- c(9.29, 12.86, 9.73, 11.45, 10.13, 9.55, 9.00, 9.78,
       12.74, 9.49, 9.70, 13.38, 9.08, 13.35, 9.33, 14.31,
       10.10, 10.03, 10.86, 9.76) #  the given dataset
p.0 <- 0.64553       # initial value
outNR <- newton(p.0, x, U, J)
outEM <- EM(x, p.0, ll)

plot(outNR$trace[,1],pch=15,col="DarkTurquoise",ylab="Estimation",main="The comparison between NR and EM algorithms")
points(outEM$trace[,1],pch=16,col="DeepPink",cex=1)
lines(outNR$trace[,1],col="DarkTurquoise",lty=1)
lines(outEM$trace[,1],col="DeepPink",lty=2)

