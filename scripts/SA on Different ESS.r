

# library(latex2exp)

# plot_dir = "Presentations/RichCon 2024/Figures/"



#
# Define the target distribution, f(x) standard Normal distribution.
# The function Takes x and returns p.d.f of N(0,1)
#
f = function(x){
  ( 1/sqrt(2*pi) ) * exp( -0.5*(x)^2 )
}


#
# Define the proposal function, g(x,\sigma) N(0, \sigma)
# A function of x (sample) and \sigma parameter (SD)
# For any values of x and \sigma, returns the p.d.f of N(0,\sigma)
#
g = function(x, sigma){
    dnorm(x, 0, sigma)
}

sim_G = function(N, sigma){
    rnorm(N, 0, sigma)
}


ESS_from_X = function(sigma, some_Xs){
    some_Ws = f(some_Xs) / g(some_Xs, sigma)
    return(1/sum(some_Ws^2))
}

ESS_from_N = function(sigma, N){
    some_Xs = sim_G(N, sigma)
    ESS_from_X(sigma, some_Xs)
}


#
# Analytical gradient of the effective sample size.
#
G = function(mu){
  2*mu*exp(mu^2)
}

# ToDo: Implement "common random number" variance reduction technique for the step-up and step-down evaluations. Switch from N as input to a set of uniforms as input. Use, e.g., qnorm(U, 0, sigma +/- c) for common random numbers
# c: step size for finite difference. Denominator is 2c
fin_diff_grad = function(sigma, N, c){
    A = ESS_from_N(sigma + c, N)
    B = ESS_from_N(sigma - c, N)

    return((A - B)/(2*c))
}

# c = get_c(2)
# fin_diff_grad(2, N, get_c(2))

#
# Define a function to update the proposal at each iteration, let's call this function update_proposal.
# This function takes the sample (x) and parameter (sigma) at each iteration and returns
# a new value for the parameter (sigma).
#
# (Some notes on stepsize:
# the best practical suggestion is to set stepsize as \frac{1}{t} where t is the iteration number.
# Usually the range for stepsize is between \frac{1}{sqrt(t)} and \frac{1}{t}
# Anything in this ranges should work it is a matter of time for convergence.)
#
# a: step size multiplier on gradient
# c: step size for finite difference. Denominator is 2c
update_proposal = function(sigma, N, a, c){
  grad_hat = fin_diff_grad(sigma, N, c)
  
  sigma_new = sigma + a * grad_hat
  return(sigma_new)
}





# # Kiefer-Wolfowitz
# get_a = function(k) return(k^(-1))
# get_c = function(k) return(0.001 * k^(-1/3))

# # Aggressive
# get_a = function(k) return(k^(-0.8))
# get_c = function(k) return(k^(-0.25))

# Chen et al. (2024)
get_a = function(k) return(k^(-0.501))
get_c = function(k) return(0.001 * k^(-0.501))



set.seed(111)

#
# Initialize values
#

# Number of SA iterations
# MC = 100
MC = 1000

# Initial value for sigma in proposal
sigma = numeric(MC)

sigma[1] = 2


# Number of random sample at each iteration
# N = 100
N = 1000

for(i in 2:MC){
  print(paste0(i, " out of ", MC))
  # You can change stepsize if you run into problem
  # but it should be ok with this one (William: I think so).
  sigma[i]  = update_proposal(sigma[i-1], N, get_a(i), get_c(i) )
}


# ToDo: Make a plot of the true ESS as a function of sigma. Compute the true value by Monte Carlo with high precision



# 
# 
# pdf(paste0(plot_dir, "PS traj.pdf"), width=10, height=7)
# 
# par(mfrow=c(1,2))
# # par(mfrow=c(1,3))
# 
# plot(1:MC, mu, xlab = 'Iteration', ylab = TeX(r'($\hat{\theta}$)'), main = 'Parameter Estimate')
# # abline(h = 0)
# plot(1:MC, all_k_hats, xlab = 'Iteration', ylab = TeX(r'($\hat{k}$)'), main = 'Tail Index')
# 
# dev.off()

# plot(20:MC, cumsum(mu[20:MC]) / 20:MC)#, ylim = c(0, 1.5))
# 
# 
# mean(mu[20:MC])

mu[MC]
mean(mu[start:MC])



MC_small = 100


pdf(paste0(plot_dir, "PS traj.pdf"), width=10, height=7)
par(mfrow=c(1,2))

plot(1:MC_small, mu[1:MC_small], xlab = 'Iteration', ylab = TeX(r'($\hat{\theta}$)'), main = 'Parameter Estimate')
# abline(h = 0)
plot(1:MC_small, all_k_hats[1:MC_small], xlab = 'Iteration', ylab = TeX(r'($\hat{k}$)'), main = 'Tail Index', ylim = c(-0.5, 0.5))

dev.off()


pdf(paste0(plot_dir, "PS mean traj.pdf"), width=5, height=7)

# par(mfrow=c(1,2))

start = MC/2

par(mfrow = c(1,1))
# plot(start:MC, mu[start:MC], xlab = 'Iteration', ylab = TeX(r'($\hat{\theta}$)'), main = 'Parameter Estimate')
plot(start:MC, cumsum(mu[start:MC])/1:(MC - start + 1), xlab = 'Iteration', ylab = "Cumulative Average", main = 'PS - Based', ylim = c(-6e-4, 6e-4))
# abline(h = 0)

dev.off()


mu[MC_small]
mu[MC]
mean(mu[start:MC])

