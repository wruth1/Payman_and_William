#
# Implementation of stochastic gradiant optimised adaptive importance sampling (normalised case)
# Description:
# This script implements the normalised case of stochastic gradiant optimised adaptive importance 
# sampling method based on the paper: 
# "Convergence rates for optimised adaptive importance samplers" by Omer Deniz Akyildiz
# and Joaquin Miguez in Statistics and Computing Journal. 
#
# Details: Algorithm 2, Stochastic gradient OAIS on page 7 of 17 of the paper. 
#
# For this example, we work on the following example in this script: 
# target:    f(x) X ~ Gamma(\alpha,\beta)
# proposal:  g(x) X ~ Exponential(\lambda)
#


using Distributions
using Statistics
using SpecialFunctions
using Plots


#
# Define the target distribution, f(x) Gamma distribution.
# The function Takes x, alpha, and beta and returns p.d.f of Gamma(\alpha,\beta)
# 
function f(x, alpha, beta)
    (beta^(alpha))/(gamma(alpha)) * x.^(alpha-1) .* exp.(-beta * x)
end


#
# Define the proposal function, g(x,lambda) Exponential(\lambda)
# A function of x (sample) and lambda parameter
# For any values of x and \lambda, it returns the p.d.f of Exponential(\lambda)
#
function g(x, lambda)
    lambda * exp.(-lambda * x)
 end


 #
 # Analytical gradient of the effective sample size.
 #
 function compute_exact_gradient()
    # Still working on it
 end


 #
 # Unbiased estimator for the gradient of effective sample size.
 #
 # Not complete, check for alpha here.
 function estimate_gradient(x, lambda)
    ( (-1/lambda) + x ) .* f(x, alpha, beta).^2 .* g(x, lambda).^2
 end



#
# Define a function to update the proposal at each iteration, let's call this function update_proposal.
# This function takes the sample (x) and parameters (\alpha and \beta) at each iteration and returns
# a new value for the parameters. 
# 
function update_proposal(x, alpha, beta, lambda, stepsize)

   num = f.(x,alpha,beta).^2
   den = g.(x,lambda).^2
   temp = -( (1/lambda) .+ x) .* (num ./ den)
   gt = mean(temp)
   theta = [alpha,beta]
   new_theta = theta .- gt*stepsize
   return new_theta

end

#
# Define a function to compute weights at each iteration.
# The inputs are sample (x) and parameter (mu) at each iteration.
# The ouput is a vector of weights for each observation in x.
#
function compute_weights(x, alpha, lambda, normalize = true)
    weights = f.(x,alpha,lambda) ./ g.(x,lambda)
    if normalize
        weights = weights ./ mean(weights)
    end
    return weights
end


#
# A function to sample from the proposal distribution, takes sample size (n) and 
# scale in exponential distribution (lambda). Returns a vector of length n.
#
function sample_from_proposal(n,lambda)
    u = rand(Exponential(lambda),n)
    return u
end



# Initial value for lambda in proposal distribution
lambda = 1
shape = 1
scale = 2

# Number of Monte Carlo samples
T = 100

# Number of random sample at each iteration
N = 1000

xx = sample_from_proposal(N,lambda)
for t in 1:T
    theta = update_proposal(xx, shape, scale, lambda, 0.2*(1/t)^(0.55))
    xx    = sample_from_proposal(N,theta[2])
    ww    = compute_weights(xx,theta[1],theta[2])
    println(theta)
end

histogram(xx)
