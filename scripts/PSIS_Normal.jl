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
using Plots
using LaTeXStrings  # For strings of the form L"...", which are interpreted using Latex syntax. Math mode is used for the whole string unless you include $'s.


# Load the functions from Pareto_Smoothing.jl
include("./src/Pareto_Smoothing.jl")




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

# List of mu values to try
mu_vec = collect(-5:0.1:5)


# Size of each sample
N = 1000

# Number of samples to generate at each value of mu
M = 500

all_k_hats = zeros(length(mu_vec), M)


using Suppressor

@suppress_err for (i, mu) in enumerate(mu_vec)
    for m in 1:M
        xx = sample_from_proposal(N,mu)
        ww = compute_weights(xx,mu, false)
        _, k_hat = pareto_smooth(ww, M_keep="default", return_k=true)
        all_k_hats[i,m] = k_hat
    end
end



# Plot the mean k_hat for each value of mu
plot(mu_vec, mean(all_k_hats, dims=2), xlabel=L"\mu", ylabel=L"\hat{k}", label=L"mean $\hat{k}$", title=L"Mean $\hat{k}$ for each value of $\mu$")
# Add a horizontal line at 0.7
hline!([0.7], label=nothing, color=:red, linestyle=:dash)

# Add confidence bands
ucls = mean(all_k_hats, dims=2) .+ 1.96*std(all_k_hats, dims=2) / sqrt(M);
lcls = mean(all_k_hats, dims=2) .- 1.96*std(all_k_hats, dims=2) / sqrt(M);

plot!(mu_vec, ucls, label="95% CI", color=:black, linestyle=:dash)
plot!(mu_vec, lcls, label=nothing, color=:black, linestyle=:dash)

