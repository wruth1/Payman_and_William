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


using Distributions
using Statistics







#
# Define the target distribution, f(x) standard Normal distribution.
# The function Takes x and returns p.d.f of N(0,1)
# 
function f(x)
    ( 1/sqrt(2*pi) ) .* exp.( -0.5*(x).^2 )
end


#
# Define the proposal function, g(x,mu) N(\mu,1) 
# A function of x (sample) and mu parameter (mean)
# For any values of x and \mu, returns the p.d.f of N(\mu,1)
#
function g(x, mu)
    (1/sqrt(2*pi)) .* exp.(-0.5 * (x.-mu).^2 )
 end


 #
 # Analytical gradient of the effective sample size.
 #
 function G(mu)
    2*mu*exp(mu^2)
 end

function get_weight(x, mu)
    return f(x) / g(x, mu)
end

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
function update_proposal(x, mu, stepsize)

   num = f.(x).^2
   den = g.(x,mu).^2
   temp = (mu .- x) .* (num ./ den)
   gt = mean(temp)
   new_mu = mu - gt*stepsize
   return new_mu

end


# Update the proposal using our analytical expression for the gradient of the effective sample size. Note that this analytical update does not require any Monte Carlo samples.
function update_proposal_analytical(mu, stepsize)
    new_mu = mu - G(mu)*stepsize
    return new_mu
end

#
# Define a function to compute weights at each iteration.
# The inputs are sample (x) and parameter (mu) at each iteration.
# The ouput is a vector of weights for each observation in x.
#
function compute_weights(x, mu, normalize=true)
    weights = f.(x) ./ g.(x,mu)
    if normalize
        weights = weights ./ mean(weights)
    end
    return weights
end


function sample_from_proposal(n,mu)
    u = rand(Normal(mu,1),n)
    return u
end


#
# Initialize values
#

# Initial value for mu in proposal
mu = 1
# let mu = 1

# Number of Monte Carlo samples
T = 100

# Number of random sample at each iteration
N = 1000

xx = sample_from_proposal(N,mu)
for t in 1:T
   # You can change stepsize if you run into problem 
   # but it should be ok with this one (William: I think so).
    mu = update_proposal(xx,mu, 0.2*(1/t)^(0.55) ) 
    xx = sample_from_proposal(N,mu)
    ww = compute_weights(xx,mu)
    println(mu)
    #println(sum(ww .* xx))
end
# end



ww = compute_weights(xx,mu, false)
var(ww)

ww_norm = compute_weights(xx,mu)
var(ww_norm)



mu = 1
for t in 1:T
    # You can change stepsize if you run into problem 
    # but it should be ok with this one (William: I think so).
     mu = update_proposal_analytical(mu,0.3*(1/t)^(0.7) )
     println(mu)
     #println(sum(ww .* xx))
 end

