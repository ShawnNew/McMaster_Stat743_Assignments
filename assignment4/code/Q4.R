# Load data
x <- c(9.21, 11.51, 12.79, 11.85, 9.97, 8.79, 9.69, 9.68, 9.19)
y <- c(7.53, 7.48, 8.08, 8.09, 10.15, 8.4, 10.88, 6.13, 7.9, 7.05, 7.48, 7.58, 8.11)
n <- length(x)
m <- length(y)
theta_hat <- mean(x)/mean(y)

# infinitesimal jackknif
b_inf <- (1/n) * sum(x-mean(x)) / mean(y) + (1/m) * sum(y-mean(y)) * theta_hat / mean(y)
v_inf <- (n^(-2) * sum((x-mean(x))^2) + m^(-2) * theta_hat^2 * sum((y-mean(y))^2)) / mean(y)^2

# regular jackknife
b_jack <- (theta_hat/m * sum(y-mean(y)) - (1/n) * sum(x-mean(x))) / mean(y)
v_jack <- (1/(n*(n-1)) * sum((x-mean(x))^2) + theta_hat^2/(m*(m-1)) * sum((y-mean(y))^2)) / mean(y)^2

# bootstrap
r <- 1000
f <- function(x, y, r){
  mean_x <- mean(sample(x,length(x),replace = TRUE))
  mean_y <- mean(sample(y, length(y), replace = TRUE))
  theta <- mean_x / mean_y
  for(i in 2:r) {
    a <- mean(sample(x,length(x),replace = TRUE)) / mean(sample(y, length(y), replace = TRUE))
    theta <- c(theta, a)
  }
  theta_star <- theta
}

theta_star <- f(x, y, r)
b_boot <- (1/r) * sum(theta_star - theta_hat)
v_boot <- 1 / (r-1) * sum(theta_star - mean(theta_star))

  