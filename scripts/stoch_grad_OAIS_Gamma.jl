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
# target:    f(x) X ~ Gamma(\alpha,\lambda)
# proposal:  g(x) X ~ Exponential(\lambda)
#


using Distributions
using Statistics
using SpecialFunctions
using Plots


#
# Define the target distribution, f(x) Gamma distribution.
# The function Takes x and returns p.d.f of Gamma(\alpha,\lambda)
# 
function f(x, alpha, lambda)
    (lambda^(alpha))/(gamma(alpha)) * x.^(alpha-1) .* exp.(-lambda * x)
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
 function compute_exact_gradient(alpha, lambda)
    (1+ (2*alpha - 1)/(lambda) ) * (-gamma(2*alpha-1)/lambda*gamma(alpha)*gamma(alpha))
 end

#
# Define a function to update the proposal at each iteration, let's call this function update_proposal.
# This function takes the sample (x) and parameters (\alpha and \lambda) at each iteration and returns
# a new value for the parameters. 
#
# (Some notes on stepsize: 
# the best practical suggestion is to set stepsize as \frac{1}{t} where t is the iteration number.
# Usually the range for stepsize is between \frac{1}{sqrt(t)} and \frac{1}{t}
# Anything in this ranges should work it is a matter of time for convergence.)
# 
function update_proposal(x, alpha, lambda, stepsize)

   num = f.(x,alpha,lambda).^2
   den = g.(x,lambda).^2
   temp = -( (1/lambda) .+ x) .* (num ./ den)
   gt = mean(temp)
   theta = [alpha,lambda]
   new_theta = theta .- gt*stepsize
   return new_theta

end


#
# Update the proposal using our analytical expression for the gradient of the effective sample size. 
# Note that this analytical update does not require any Monte Carlo samples.
# 
function update_proposal_analytical(theta, stepsize)
    new_theta = theta .- compute_exact_gradient(theta[1],theta[2])*stepsize
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


function sample_from_proposal(n,lambda)
    u = rand(Exponential(lambda),n)
    return u
end


#
# Initialize values
#

# Initial value for alpha and lambda in proposal distribution
alpha = 1
lambda = 1

# Number of Monte Carlo samples
T = 100

# Number of random sample at each iteration
N = 1000

xx = sample_from_proposal(N,lambda)
for t in 1:T
    theta = update_proposal_analytical( [alpha,lambda], 0.2*(1/t)^(0.55) ) 
    xx    = sample_from_proposal(N,theta[2])
    ww    = compute_weights(xx,theta[1],theta[2])
    println(theta)
    println(mean(ww .* xx))
end

histogram(xx)
