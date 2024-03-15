using Distributions
using Statistics
using SpecialFunctions

#
# Note to Payman: For every value of theta_j, you need to loop over your sample size 
# and compute the value of k. Now, you can compute l(theta) (profile log-likelihood)
# and then compute likelihood. Then repeat this for all values of j = 1,...,m
# Think about k as k_j since it is in fact changes based on the value of j
#


#
# Estimate parameters in Generalized Pareto Distribution (GDP) using the method described in:
# Jin Zhang and Michael A. Stephens (2009), section 4: new estimators
#

function fit_generalized_pareto_dist(x)
    
    # Get the sample size
    n = length(x)

    # Calculate m
    m = 20 + floor( sqrt(n) )

    # Estimate theta
    theta_hat = estimate_theta(x, m)

    # Estimate k_hat
    k_hat = - mean( log.( 1 .- theta_hat .* x) )

    # Estimate sigma_hat
    sigma_hat = k_hat / theta_hat

    return sigma_hat, k_hat, theta_hat
end

function estimate_theta(x, m)

    # Calculate first quartile of the sample
    xstar = quantile(x, 0.25)

    # Compute theta_j defined in the paper
    all_thetas = (1/maximum(x)) .+ [ 1 - sqrt(m/(j-0.5))  for j in 1:m] ./ (3 * xstar)

    # Compute profile log likelihood for all_thetas as defined in the paper
    ll_theta = compute_profile_log_likelihood(x, all_thetas)

    # Compute weights for each theta
    theta_weights = compute_theta_weights(ll_theta)

    # Compute estimat eof theta for GPD
    theta_hat = sum( ll_theta .* theta_weights)
 
    return theta_hat
end

function compute_profile_log_likelihood(x, theta)

   # Get the sample size
   n = length(x)
   ll_theta= zeros(length(theta))

   for i in eachindex(theta)
      # Compute k_j based on theta_j to be able to compute profile likelihood
      k = - mean( log.(1 .- x .* theta[i]) )

     # Compute profile log-likelihood
     ll_theta[i] = n * ( log(theta[i]/k) + k - 1 )
   end

   return ll_theta
end

function compute_theta_weights(ll_theta)
       
    sm = sum( exp.(ll_theta) )
    theta_weights = (1/sm) .* ( 1 ./ exp.(-ll_theta) )
    return theta_weights

end 

x = rand(Exponential(3), 30)

fit_generalized_pareto_dist(x)

fit_GPD(x)