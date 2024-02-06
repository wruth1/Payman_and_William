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
 function compute_exact_gradient(lambda, alpha, beta)
    # Still working on it
    num = 2 * beta^(2*alpha) * gamma(2*alpha - 1) * (lambda - beta)
    den = lambda^2 * gamma(alpha)^2 * (2*alpha)
 end


 #
 # Unbiased estimator for the gradient of effective sample size.
 #
 # Not complete, check for alpha here.
 function estimate_gradient(x, lambda, alpha, beta)
    ( (-1/lambda) + x ) .* f(x, alpha, beta).^2 ./ g(x, lambda).^2
 end


 # Ref(parameter) tells Julia to replace the vector to scaler point to it. broadcast of a vector in a function

#
# Define a function to update the proposal at each iteration, let's call this function update_proposal.
# This function takes the sample (x) and parameters (\alpha and \beta) at each iteration and returns
# a new value for the parameters. 
# 
function update_proposal(x, alpha, beta, lambda, stepsize)

   temp = estimate_gradient.(x, lambda, alpha, beta)
   gt = mean(temp)
   new_lambda = lambda .- gt*stepsize
   return new_lambda

end

#
# Define a function to compute weights at each iteration.
# The inputs are sample (x) and parameter (mu) at each iteration.
# The ouput is a vector of weights for each observation in x.
#
function compute_weights(x, lambda, alpha, beta, normalize = true)
    weights = f.(x,alpha,beta) ./ g.(x,lambda)
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

# Fix the value of shape and scale in Gamma dist
shape = 1
scale = 5


# Initial value for lambda in proposal distribution
lambda_0 = 4.0
lambda = lambda_0
theta = lambda_0

# Number of Monte Carlo samples
T = 5000

# Number of random sample at each iteration
N = 2000

all_thetas = zeros(T)
all_ess = zeros(T)

xx = sample_from_proposal(N,lambda)
for t in 1:T
    theta = update_proposal(xx, shape, scale, theta, 0.5*(1/t)^(0.9))
    xx    = sample_from_proposal(N,theta)
    ww    = compute_weights(xx,theta,shape,scale)
    all_thetas[t] = theta
    all_ess[t] = mean( f.(xx,shape,scale).^2 / g.(xx,scale) )
    println(theta)
end

scatter(all_thetas)
scatter(all_ess)
#ylims!(lambda_0*0.95, maximum(all_thetas))
#hline!([lambda_0])
