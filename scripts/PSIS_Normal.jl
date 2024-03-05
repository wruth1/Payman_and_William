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
 # Analytical gradient of the effective sample size.
 #
 function G(mu)
    2*mu*exp(mu^2)
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


"""
Analytical value of the objective function we are trying to minimize.
"""
function rho(mu)
    exp(mu^2)
end

"""
Estimate of our objective function obtained from an importance sample.
"""
function rho_hat(all_Xs, mu)
    all_wts = compute_weights(all_Xs, mu, false)
    return mean(all_wts .^2)
end


#
# Initialize values
#

# List of mu values to try
## Note: rho grows like exp(mu^2), so you really have to zoom-in to see what's happening near zero
# mu_vec = collect(-3:0.1:3)
mu_vec = collect(-2:0.1:2)


# Size of each sample
N = 1000

# Number of samples to generate at each value of mu
M = 500

all_k_hats = zeros(length(mu_vec), M)
all_rho_hats = zeros(length(mu_vec), M)

using Suppressor

@suppress_err for (i, mu) in enumerate(mu_vec)
    for m in 1:M
        xx = sample_from_proposal(N,mu)
        ww = compute_weights(xx,mu, false)
        _, k_hat = pareto_smooth(ww, M_keep="default", return_k=true)
        all_k_hats[i,m] = k_hat

        all_rho_hats[i,m] = rho_hat(xx, mu)
    end
end

all_true_rhos = rho.(mu_vec)



# Plot the mean k_hat for each value of mu
## Use a hacky solution to use a single legend. I will add this label later.
plot(mu_vec, mean(all_k_hats, dims=2), xlabel=L"\mu", ylabel=L"\hat{k}", label=nothing, title=L"Mean $\hat{k}$ for each value of $\mu$")
# Add a horizontal line at 0.7
hline!([0.7], label=nothing, color=:red, linestyle=:dash)

# Add confidence bands
ucls = mean(all_k_hats, dims=2) .+ 1.96*std(all_k_hats, dims=2) / sqrt(M);
lcls = mean(all_k_hats, dims=2) .- 1.96*std(all_k_hats, dims=2) / sqrt(M);

#* Note: I remove the confidence bands because having two sets of them makes the plot too messy.
# plot!(mu_vec, ucls, label=nothing, color=:black, linestyle=:dashdot)
# plot!(mu_vec, lcls, label=nothing, color=:black, linestyle=:dashdot)



# Add estimated and true values of rho. Use the right vertical axis.
rhs = twinx()
plot!(rhs, mu_vec, all_true_rhos, ylabel=L"\rho(\mu)", label=L"\rho", color=:green, linestyle=:dash)

# Add mean estimated values of rho. Use the right vertical axis.
plot!(rhs, mu_vec, mean(all_rho_hats, dims=2), label=L"\hat{\rho}", color=:purple, linestyle=:dash)

# Add confidence bands
ucls_rho = mean(all_rho_hats, dims=2) .+ 1.96*std(all_rho_hats, dims=2) / sqrt(M);
lcls_rho = mean(all_rho_hats, dims=2) .- 1.96*std(all_rho_hats, dims=2) / sqrt(M);

#* Note: I remove the confidence bands because having two sets of them makes the plot too messy.
# plot!(rhs, mu_vec, ucls_rho, label=nothing, color=:black, linestyle=:dashdot)
# plot!(rhs, mu_vec, lcls_rho, label=nothing, color=:black, linestyle=:dashdot)


# Add legend for original plot
plot!(rhs, mu_vec, NaN.*mean(all_k_hats, dims=2), xlabel=L"\mu", ylabel=L"\hat{k}", label=L"mean $\hat{k}$", title=L"Mean $\hat{k}$ for each value of $\mu$")


