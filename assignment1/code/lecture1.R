# Stat 743B Winter 2018

# Code for examples done in class on January 4 and 8, 2018

# Generating an exponential with rate equal to 10 using the
# inverse method
u1 <- runif(100000)
y1 <- -log(1-u)/10
hist(u1, breaks=100, prob=T)
abline(h=1, col="red")
hist(y1, breaks=100, prob=T)
lines(seq(0,max(y1), by=0.01), dexp(seq(0, max(y1), by=0.01),10),
      col="red", lwd=2)

# Generating a Binomial(n=3, p=0.5) using the inverse method.
u2 <- runif(100000)
y2 <- rep(0,length(u))
y2[u2>0.125&u2<=0.500] <- 1
y2[u2>0.5&u2<=0.875] <- 2
y2[u2>0.875] <- 3
table(y2)/length(y2)

# Accept-Reject example to generate a beta(1.5,3) using a uniform candidate
M <- dbeta(0.2,1.5,3)
V <- runif(100000)
U <- runif(100000)
prob <- dbeta(V,1.5,3)/(M*dunif(V))
Accept <- U<prob
Y <- V[Accept]
length(Y)
hist(Y, breaks=100, prob=T)
lines((0:1000)/1000, dbeta((0:1000)/1000, 1.5,3), col="red", lwd=2)