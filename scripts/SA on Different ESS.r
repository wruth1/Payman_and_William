
# ---------------------------------------------------------------------------- #
#                         Easy - Normal Location Family                        #
# ---------------------------------------------------------------------------- #

set.seed(1)

source("src/ESS Versions.r")

N = 1000

mu0 = 1

X = rnorm(N, mu0, 1)


mu = mu0


# Kiefer-Wolfowitz
get_a = function(k) return(k^(-1))
get_c = function(k) return(k^(-1/3))

# # Aggressive
# get_a = function(k) return(k^(-0.8))
# get_c = function(k) return(k^(-0.25))

wt = function(theta, X, normalize = T){
    W_raw = dnorm(X) / dnorm(X, theta, 1)

    if(!normalize){
        return(W_raw)
    } else{
        W = W_raw / sum(W_raw)
        return(W)
    }    
} 

ESS = function(theta, X){
    W = wt(theta, X)

    # L2
    return(1 / sum(W^2))

    # L Infinity
    # return(1 / max(W))
    
}

rproposal = function(n, theta) rnorm(n, theta, 1)


fin_diff_ESS = function(theta, X, c){
    A = ESS(theta + c, X)
    B = ESS(theta - c, X)

    return((A - B) / (2*c))
}



# k: iteration number
# n: size of sample from proposal
update_theta = function(theta_old, k, n=1000){
    a_k = get_a(k)
    
    X = rproposal(n, theta_old)

    c_k = get_c(k)
    grad_hat = fin_diff_ESS(theta_old, X, c_k)

    theta_new = theta_old + a_k * grad_hat
    return(theta_new)
}



set.seed(1)

K = 10

theta_old = 1

theta_traj = rep(0, times = K)
theta_traj[1] = theta_old

for(k in 2:K){
    theta_next = update_theta(theta_old, k, n)

    theta_traj[k] = theta_next
    theta_old = theta_next
}

theta_traj
