
# set.seed(1)

# # Target: standard normal
# # Proposal: N(1,2)

# mu_prop = 1
# sd_prop = 2

# n = 100

# X = rnorm(n, mu_prop, sd_prop)

# W = dnorm(X) / dnorm(X, mu_prop, sd_prop)

# W_norm = W / sum(W)


# Effective sample size based on L2 distance from discrete unif.
# Assumes weights are already self-normalized
ESS_L2 = function(W){
    return(1 / sum(W^2))
}

ESS_LInf = function(W){
    1/max(W)
}

ESS_L1 = function(W){
    n = length(W)
    Wp = subset(W, W > 1/n)

    Np = length(Wp)
    Sp = sum(Wp)

    return(-n * Sp + Np + n)
}

# ESS_L2(W_norm)
# ESS_LInf(W_norm)
# ESS_L1(W_norm)


