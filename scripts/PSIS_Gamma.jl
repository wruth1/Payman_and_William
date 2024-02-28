
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


using Distributions
using Statistics
using Plots
using LaTeXStrings  # For strings of the form L"...", which are interpreted using Latex syntax. Math mode is used for the whole string unless you include $'s.


# Load the functions from Pareto_Smoothing.jl
include("./src/Pareto_Smoothing.jl")






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




#
# A function to sample from the proposal distribution, takes sample size (n) and 
# scale in exponential distribution (lambda). Returns a vector of length n from exponential distribution.
#
function sample_from_proposal(n,lambda)
    u = rand(Exponential(lambda),n)
    return u
end


#
# Initialize values
#

# Fix the value of shape and scale in Gamma distribution (our target distribution)
shape = 1
scale = 2


# Initial value for lambda in proposal distribution 
lambda = 1
lambda_vec = 0.5:0.1:5


# Size of each sample
N = 1000

# Number of samples to generate at each value of mu
M = 50

all_k_hats = zeros(length(lambda_vec), M)


using Suppressor

@suppress_err for (i, lambda) in enumerate(lambda_vec)
    for m in 1:M
        xx = sample_from_proposal(N,lambda)
        ww = compute_weights(xx,lambda, false)
        _, k_hat = pareto_smooth(ww, M_keep="default", return_k=true)
        all_k_hats[i,m] = k_hat
    end
end



# Plot the mean k_hat for each value of lambda
plot(lambda_vec, mean(all_k_hats, dims=2), xlabel=L"\lambda", ylabel=L"\hat{k}", label=L"mean $\hat{k}$", title=L"Mean $\hat{k}$ for each value of $\lambda$")
# Add a horizontal line at 0.7
hline!([0.7], label=nothing, color=:red, linestyle=:dash)

# Add confidence bands
ucls = mean(all_k_hats, dims=2) .+ 1.96*std(all_k_hats, dims=2) / sqrt(M);
lcls = mean(all_k_hats, dims=2) .- 1.96*std(all_k_hats, dims=2) / sqrt(M);

plot!(lambda_vec, ucls, label="95% CI", color=:black, linestyle=:dash)
plot!(lambda_vec, lcls, label=nothing, color=:black, linestyle=:dash)

#! This is weird. Negative values of k-hat are possible (if I remember correctly), but I'm not really sure how to interpret them beyond "something is probably wrong". In particular, I don't think Vehtari et al. every discuss negative k-hats.
#! I'm going to add the ESS-type objective function on the same plot so we can compare. I'm not sure if we ever plotted that objective vs lambda. If not, then this should be instructive.