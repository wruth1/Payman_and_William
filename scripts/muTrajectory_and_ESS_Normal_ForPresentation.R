#
# Implementation of stochastic gradiant optimised adaptive importance sampling (normalised case)
# Description:
# This script implements the normalised case of  stochastic gradiant optimised adaptive importance
# sampling method based on the paper:
# "Convergence rates for optimised adaptive importance samplers" by Omer Deniz Akyildiz
# and Joaquin Miguez in Statistics and Computing Journal.
#
# Details: Algorithm 2, Stochastic gradient OAIS on page 7 of 17 of the paper.
#
# I used the following example in this script:
# target:    f(x) = ( 1/sqrt(2*pi) ) exp( -0.5 x^2 )        X ~ N(0,1)
# proposal:  g(x) = ( 1/sqrt(2*pi) ) exp( -0.5 (x-\mu)^2 )  X ~ N(\mu,1)
#



#
# Define the target distribution, f(x) standard Normal distribution.
# The function Takes x and returns p.d.f of N(0,1)
#
f = function(x){
  ( 1/sqrt(2*pi) ) * exp( -0.5*(x)^2 )
}


#
# Define the proposal function, g(x,mu) N(\mu,1)
# A function of x (sample) and mu parameter (mean)
# For any values of x and \mu, returns the p.d.f of N(\mu,1)
#
g = function(x, mu){
  (1/sqrt(2*pi)) * exp(-0.5 * (x-mu)^2 )
}


#
# Analytical gradient of the effective sample size.
#
G = function(mu){
  2*mu*exp(mu^2)
}

#
# Define a function to update the proposal at each iteration, let's call this function update_proposal.
# This function takes the sample (x) and parameter (mu) at each iteration and returns
# a new value for the parameter (mu).
#
# (Some notes on stepsize:
# the best practical suggestion is to set stepsize as \frac{1}{t} where t is the iteration number.
# Usually the range for stepsize is between \frac{1}{sqrt(t)} and \frac{1}{t}
# Anything in this ranges should work it is a matter of time for convergence.)
#
update_proposal = function(x, mu, stepsize){
  num = f(x)^2
  den = g(x,mu)^2
  temp = (mu - x) * (num / den)
  gt = mean(temp)
  new_mu = mu - gt*stepsize
  return(new_mu)
}


# Update the proposal using our analytical expression for the gradient of the effective sample size. Note that this analytical update does not require any Monte Carlo samples.
update_proposal_analytical = function(mu, stepsize){
  new_mu = mu - G(mu)*stepsize
  return(new_mu)
}

#
# Define a function to compute weights at each iteration.
# The inputs are sample (x) and parameter (mu) at each iteration.
# The ouput is a vector of weights for each observation in x.
#
compute_weights = function(x, mu){
  weights = f(x) / g(x,mu)
  return(weights)
}


sample_from_proposal = function(n, mu){
  u = rnorm(n = n, mean = mu, sd = 1)
  return(u)
}


set.seed(1)

#
# Initialize values
#

# Number of Monte Carlo samples
MC = 100

# Initial value for mu in proposal
mu = numeric(MC)
ess = numeric(MC)

mu[1] = 1.5

# let mu = 1

# Number of random sample at each iteration
N = 100

xx = sample_from_proposal(N, mu[1])

for(i in 2:MC){
  # You can change stepsize if you run into problem
  # but it should be ok with this one (William: I think so).
  mu[i]  = update_proposal(xx, mu[i-1], 0.2*(1/i)^(0.55) )
  xx     = sample_from_proposal(N, mu[i])
  temp   = mean( f(xx) / g(xx, mu[i]) )^2
  ess[i] = mean( f(xx)^2 / g(xx, mu[i])^2 )  / temp
}

par(mfrow=c(1,2))
plot(x = 1:MC, y = mu, xlab = 'Iteration', ylab = 'mu', main = 'Trajectory of sampling for the mean parameter')
plot(x = 1:MC, y = N/ess, xlab = 'Iteration', ylab = 'ESS', main = 'Effective sample size')


