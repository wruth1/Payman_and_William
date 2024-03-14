using Distributions
using Statistics


#
# Estimate parameters in Generalized Pareto Distribution (GDP) using the method described in:
# Jin Zhang and Michael A. Stephens (2009), section 4: new estimators
#

function fit_gpd(x)
    
    # Get the sample size
    n = length(x)

    # Calculate m
    m = 20 + floor( sqrt(n) )

    # Estimate theta
    theta_hat = estimate_theta(x, m)
end


function estimate_theta(x, m)

    # Calculate first quartile of the sample
    xstar = quantile(x, 0.25)

    # Compute theta_j defined in the paper
    all_thetas = (1/maximum(x)) .+ [ 1 - sqrt(m/(j-0.5))  for j in 1:m] ./ (3 * xstar)

    # Compute profile likelihood for all_thetas as defined in the paper
    tmp = compute_profile_likelihood(x[1:trunc(Int,m)], all_thetas)
    pl_all_thetas = tmp / sum(tmp)

    # Compute the estimate of theta in GPD
    theta_hat = sum(all_thetas * pl_all_thetas)

    return theta_hat
end


function compute_profile_likelihood(x, theta)

   # Get the sample size
   n = length(x)

   # Compute k for the profile likelihood formula
   k = (-1/n) * ( log(1- theta'x) )

   # Compute profile log-likelihood
   l_theta = n .* ( log.(theta/k) .+ k .- 1 )

   # Compute profile likelihood and return
   L_theta = exp.(l_theta)

   return L_theta
end


x = rand(Exponential(3), 30)
fit_gpd(x)