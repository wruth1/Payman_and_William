M <- 1000

# Target distribution: N(0, 1)
p = function(x){
  return(dnorm(x))
}

# Proposal distribution: N(mu, 1)
q = function(x, mu){
  return(dnorm(x, mean = mu, sd = 1))
}

set.seed(123)
par(mfrow=c(1,2))

# Proposal 1: N(0.1, 1)
x1 <- rnorm(n = M, mean = 0.1, sd = 1)
w1 <- p(x1) / q(x1, mu = 0.1)
#hist(w1, breaks = 10, main = "Weights for Proposal N(0.1,1)", xlab = "Weights")
hist(w1, 'fd', main = 'Weights for Proposal N(mean=0.1,sd=1)', xlab = 'Weights')


# Proposal 2: N(2, 1)
x2 <- rnorm(n = M, mean = 2, sd = 1)
w2 <- p(x2) / q(x2, mu = 2)
hist(w2, breaks = 100, main = 'Weights for Proposal N(mean=2,sd=1)', xlab = 'Weights')


