# R code to implement the Newton-Raphson algorithm
#


U <- function(pi, x)
  # The data x are ten-dimensional
  # The derivative of the log-likelihood for pi 
  sum((dnorm(x[1, ], mean=10, sd=1) - dnorm(x[2, ], mean=13, sd=1)) / 
        (pi * dnorm(x[1, ], mean=10, sd=1) + (1-pi) * dnorm(x[2, ], mean=13, sd=1)))
  #log(al)-digamma(al)-log(mean(x))+mean(log(x))

J <- function(pi, x) 
  # The data x are ten-dimensional
  # The NEGATIVE of the second derivative of the log-likelihood
  # for pi
  -sum((dnorm(x[, 1], mean=10, sd=1) - dnorm(x[, 2], mean=13, sd=1))^2 / 
         (pi * dnorm(x[, 1], mean=10, sd=1) + (1-pi) * dnorm(x[, 2], mean=13, sd=1))^2)
  #trigamma(al)-1/al

newton <- function (th0, x, U, J, eps=1e-12, maxit=1000000) {
  # A general function to implement a ten-dimensional 
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
    continue <- ((abs(out[iter,2])>eps) | 
                   (out[iter,4])>eps) & (iter<=maxit))
  }
  out <- out[1:iter,]
  return(list(est=out[iter,1], trace=out))
}

# Here is a simulated example
x <- c(9.29, 12.86, 9.73, 11.45, 10.13, 9.55, 9.00, 9.78,
       12.74, 9.49, 9.70, 13.38, 9.08, 13.35, 9.33, 14.31,
       10.10, 10.03, 10.86, 9.76) #  the given dataset
x <- matrix(x, nrow=2, ncol=10)
#alpha.0 <- mean(x)^2/var(x)
pi.0 <- 1       # initial value
output <- newton(pi.0, x, U, J)
pi.hat <- newton(pi.0, x, U, J)$est
#beta.hat <- mean(x)/alpha.hat
