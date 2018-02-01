library("MASS")
library("fMultivar")
library("timeDate")
library("timeSeries")
library("fBasics")
library("mixtools")
library("plot3D")

# Parameters for bivariate normal distribution
mu1 = 0; mu2 = 0
s1 = 1; s2 = 1
#mu <- c(0,0)
#Sigma <- matrix(c(1,0.7,0.7,1), 2,2)
rho = 0.7

# Initialization
Nsim = 10^4


# Function to draw ellipse for bivariate normal data
draw3DPic <- function(bvn, bins) {
  ##  Create cuts:
  x_c <- cut(bvn[,1], bins)
  y_c <- cut(bvn[,2], bins)
  
  ##  Calculate joint counts at cut levels:
  z <- table(x_c, y_c)
  
  ##  Plot as a 3D histogram:
  hist3D(z=z, border="black")
}



# Function of Independent Metropolis-Hastings algorithm
iMH <- function (n, mu1, s1, mu2, s2, rho) {
  # Parameters for bivariate normal distribution
  mu <- c(mu1,mu2) # Mean 
  Sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2, 2) # Covariance matrix
  I <- matrix(c(1,0,0,1),2,2)
  mat <- matrix(nrow = n, ncol = 2)
  mat[1, ] <- mvrnorm(1, mu, Sigma)    # initialize the chain from the stationary
  for (t in 2 : n) {
    # Sample from proposal
    xStar <- mvrnorm(1, mat[t-1, ], I)
    
    # The acceptance ratio
    c = dmvnorm(mat[t-1, ], xStar) / dmvnorm(xStar, mat[t-1, ])
    alpha = min(1, dmvnorm(xStar) / dmvnorm(mat[t-1, ]))
    
    # Decide to accept
    mat[t, ] <- mat[t-1, ] + (xStar - mat[t-1, ]) * (runif(1)<alpha)
  }
  mat
}

# Function of random walk Metropolis-Hastings algorithm








# Function of Gibbs Sampler to generate
gibbs<-function (n, mu1, s1, mu2, s2, rho) 
{
  mat <- matrix(ncol = 2, nrow = n)
  x <- 0
  y <- 0
  mat[1, ] <- c(x, y)
  for (i in 2:n) {
    x <- rnorm(1, mu1 + 
                 (s1/s2) * rho * (y - mu2), sqrt((1 - rho^2)*s1^2))
    y <- rnorm(1, mu2 + 
                 (s2/s1) * rho * (x - mu1), sqrt((1 - rho^2)*s2^2))
    mat[i, ] <- c(x, y)
  }
  mat
}