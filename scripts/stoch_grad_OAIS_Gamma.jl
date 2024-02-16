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
# For this example, we work on the following target and proposal distribution in this script: 
# target:    f(x) X ~ Gamma(\a,\beta)
# proposal:  g(x) X ~ Exponential(\lambda)
#


#
# Load the required packages
#
using Distributions
using Statistics
using SpecialFunctions
using Plots
using Random
using StatsPlots # to save plots
using ProgressMeter



#
# Define the target distribution, f(x,alpha,beta) Gamma distribution.
# The function takes x, alpha (shape parameter), and beta (scale parameter) 
# and returns p.d.f of Gamma(\alpha,\beta)
# 
function f(x, alpha, beta)
    res = (beta^(alpha))/(gamma(alpha)) * x.^(alpha-1) .* exp.(-beta * x)
    return(res)
end


#
# Define the proposal distribution, g(x,lambda) Exponential(\lambda)
# The function takes x and lambda (scale parameter in exp distribution) and 
# returns p.d.f of Exponential(\lambda)
#
function g(x, lambda)
    res = lambda * exp.(-lambda * x)
    return(res)
 end


 #
 # Analytical calculation of gradient of the effective sample size.
 #
 function compute_exact_gradient(lambda, alpha, beta)

    term1 = ( (beta)/(2*beta - lambda) )^(2*alpha)
    term2 = 1 / (lambda * gamma(alpha)^2 )
    term3 = gamma(2*alpha) - ( ( (2*beta - lambda)*(gamma(2*alpha - 1)) )/lambda )
    res = term1 * term2 * term3

    return(res)

 end


 #
 # Unbiased estimator for the gradient of effective sample size.
 #
 function estimate_gradient(x, lambda, alpha, beta)
    temp = f(x, alpha, beta).^2 ./ g(x, lambda).^2
    res = ( x .- (1/lambda) ) .* temp
    return(res)
 end



#
# Define a function to update the proposal at each iteration, let's call this function update_proposal.
# This function takes the followings:
#  x: sample 
# \a: shape in target distribution (constant) 
# \b: scale in target distribution (constant) 
# stepsize: the step size in optimization algorithm. new_parameter = old_parameter - (stepsize * gradient)
# 
# And returns:
# 
# new_lambda: a scalar value. The updated value for parameter in proposal distribution
# 
function update_proposal(lambda, stepsize, gt)

   new_lambda = lambda .- (stepsize * gt)
   return new_lambda

end


#
# Define a function to compute gradient either exactly or empirically
#
#  x: sample 
# \alpha: shape in target distribution (constant) 
# \beta: scale in target distribution (constant)
# exact: a boolean value, true (default value) computes the exact gradient. 
#
function compute_gradient(x, alpha, beta, lambda, exact=true)

    if exact
        # Compute the exact gradient analytically 
        gt = compute_exact_gradient(lambda, alpha, beta)
    else
        # Estimate the gradient empirically 
        temp = estimate_gradient.(x, lambda, alpha, beta)
        gt = mean(temp)
    end

    return gt

end


#
# Define a function to compute weights at each iteration.
# The functions takes the followins:
#
# x: sample
# lambda: scale parameter in exponential distribution
# a: shape value in target distribution
# b: scale value in proposal distribution
# normalize: a binary value to indicate whether we need to normalize the weights. 
# 
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
# scale in exponential distribution (lambda). Returns a vector of length n from exponential distribution.
#
function sample_from_proposal(n,lambda)
    u = rand(Exponential(lambda),n)
    return u
end


# Fix the value of shape and scale in Gamma distribution (our target distribution)
shape = 1
scale = 2


# Initial value for lambda in proposal distribution (We use theta here to be consistent with paper notation)
theta = 1

# Number of Monte Carlo samples
T = 1000

# Number of random sample at each iteration
N = 10000


#
# Define two numeric vector to record: 1) value for lambda, and 2) effective sample size (ess) 
# at each iteration (1:T)
#
all_thetas_true = zeros(T)
all_thetas_approx = zeros(T)
all_ess_true = zeros(T)
all_ess_approx = zeros(T)
all_true_gradient = zeros(T)
all_approx_gradient = zeros(T)

# Set seed
Random.seed!(123)

# Sample from proposal
#xx = sample_from_proposal(N,theta)

# Iterate over 1:T, at each iteration: 
# (i) update proposal parameter, (ii) sample from proposal with the new theta value, and 
# (iii) compute weights
# Along the loop, record each generated value for theta and effective sample size.

exact=true

@showprogress for t in 1:T

    xx    = sample_from_proposal(N,theta)

    gradient = compute_gradient(xx, shape, scale, theta, exact)
    theta = update_proposal(theta, 0.5*(1/t)^(0.9), gradient)
    ess = mean( f.(xx,shape,scale).^2 / g.(xx,scale) )

    if exact
        all_true_gradient[t] = gradient
        all_thetas_true[t] = theta
        all_ess_true[t] = ess
    else
        all_approx_gradient[t] = gradient
        all_thetas_approx[t] = theta
        all_ess_approx[t] = ess
    end

    #ww    = compute_weights(xx,theta,shape,scale)

    #println(theta)
end




exact = false

@showprogress for t in 1:T

    xx    = sample_from_proposal(N,theta)

    gradient = compute_gradient(xx, shape, scale, theta, exact)
    theta = update_proposal(theta, 0.5*(1/t)^(0.9), gradient)
    ess = mean( f.(xx,shape,scale).^2 / g.(xx,scale) )

    if exact
        all_true_gradient[t] = gradient
        all_thetas_true[t] = theta
        all_ess_true[t] = ess
    else
        all_approx_gradient[t] = gradient
        all_thetas_approx[t] = theta
        all_ess_approx[t] = ess
    end

    #ww    = compute_weights(xx,theta,shape,scale)

    #println(theta)
end



#scatter(all_thetas)
#scatter(all_ess)
#ylims!(lambda_0*0.95, maximum(all_thetas))
#hline!([lambda_0])


scatter(all_thetas_true)
hline!([2.0])
scatter(all_thetas_approx)
hline!([2.0])

scatter(all_true_gradient)
scatter(all_approx_gradient)
scatter(all_true_gradient, all_approx_gradient)


true_thetas_plot = scatter(all_thetas_true);
approx_thetas_plot = scatter(all_thetas_approx);
plot(true_thetas_plot, approx_thetas_plot, layout = (1, 2), size = (1200, 600))
savefig("trueThetas_vs_approxThetas.png")

true_grad_plot = scatter(all_true_gradient);
ylims!((-0.1, 0.05))
approx_grad_plot = scatter(all_approx_gradient, markersize=0.5);
ylims!((-0.1, 0.05))
plot(true_grad_plot, approx_grad_plot, layout = (1, 2), size = (1200, 600))
savefig("trueGradients_vs_approxGradients.png")


