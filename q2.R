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
rho = 0.7

# Sampling rate
Nsim = 10^6


# Function to 3D histogram of the bivariate normal points
draw3DPic <- function(bvn, bins) {
  # Create cuts:
  x_c <- cut(bvn[,1], bins)
  y_c <- cut(bvn[,2], bins)
  
  # Calculate joint counts at cut levels:
  z <- table(x_c, y_c)
  
  # Plot as a 3D histogram:
  hist3D(z=z, breaks=100)
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
rwMH <- function (n, mu1, s1, mu2, s2, rho) {
  # Parameters for bivariate normal distribution
  mu <- c(mu1,mu2) # Mean 
  Sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2, 2) # Covariance matrix
  I <- matrix(c(1,0,0,1),2,2)
  mat <- matrix(nrow = n, ncol = 2)
  mat[1, ] <- mvrnorm(1, mu, Sigma)    # initialize the chain from the stationary
  for (t in 2:n) {
    # Sample from proposal distribution
    xStar <- mvrnorm(1, mat[t-1, ], I)
    # Do the random walk
    #Y <- dmvnorm(xStar)
    Y <- mvrnorm(1, xStar, I)
    # Compute the probobility
    prob <- min(1, dmvnorm(Y)/dmvnorm(mat[t-1, ]))
    # Decide to accept
    mat[t, ] <- mat[t-1, ] + (Y - mat[t-1, ]) * (runif(1) < prob)
  }
  mat
}


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


# Define the main function to invoke the three methods to generate data
# and also store the data into file
myMainFunction <- function() {
  bin <- 80
  print("Generate data from Independence MH algorithm...")
  ibvn <- iMH(Nsim, mu1, s1,mu2, s2, rho)
  print("Generate data from random walk MH algorithm...")
  rwbvn <- rwMH(Nsim, mu1, s1, mu2, s2, rho)
  print("Generate data from Gibbs Sampler algorithm...")
  gibbsbvn <- gibbs(Nsim, mu1, s1, mu2, s2, rho)
  df <- data.frame(ibvn, rwbvn, gibbsbvn)    # create data frame to store in csv file
  names(df) <- c("DatafromIMH", "DatafromRWMH", "DatafromGibbs")
  print("Write the data into csv file...")
  write.csv(df,"dummyData.csv")
  print("Draw pictures")
  opar <- par(no.readonly = TRUE)
  par(mfrow=c(3,1))
  
  # iMH
  par(family="mono")
  draw3DPic(ibvn, bin)
  title(main="iMH generated data",
        sub="By Chenxiao Niu")
  
  # rwMH
  par(family="mono")
  draw3DPic(rwbvn, bin)
  title(main="Random Walk generated data",
        sub="By Chenxiao Niu")
  
  # Gibbs
  par(family="mono")
  draw3DPic(gibbsbvn, bin)
  title(main="Gibbs generated data",
        sub="By Chenxiao Niu")
  
  par(opar)
}
