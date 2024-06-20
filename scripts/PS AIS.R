

library(latex2exp)

plot_dir = "Presentations/RichCon 2024/Figures/"



#
# Define the target distribution, f(x) standard Normal distribution.
# The function Takes x and returns p.d.f of N(0,1)
#
f = function(x){
  ( 1/sqrt(2*pi) ) * exp( -0.5*(x)^2 )
}


#
# Define the proposal function, g(x,mu) N(\mu,1)
# A function of x (sample) and mu parameter (mean)
# For any values of x and \mu, returns the p.d.f of N(\mu,1)
#
g = function(x, mu){
  (1/sqrt(2*pi)) * exp(-0.5 * (x-mu)^2 )
}


#
# Analytical gradient of the effective sample size.
#
G = function(mu){
  2*mu*exp(mu^2)
}

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
update_proposal = function(x, mu, stepsize){
  num = f(x)^2
  den = g(x,mu)^2
  temp = (mu - x) * (num / den)
  gt = mean(temp)
  new_mu = mu - gt*stepsize
  return(new_mu)
}


# Update the proposal using our analytical expression for the gradient of the effective sample size. Note that this analytical update does not require any Monte Carlo samples.
update_proposal_analytical = function(mu, stepsize){
  new_mu = mu - G(mu)*stepsize
  return(new_mu)
}

#
# Define a function to compute weights at each iteration.
# The inputs are sample (x) and parameter (mu) at each iteration.
# The ouput is a vector of weights for each observation in x.
#
compute_weights = function(x, mu){
  weights = f(x) / g(x,mu)
  return(weights)
}


sample_from_proposal = function(n, mu){
  u = rnorm(n = n, mean = mu, sd = 1)
  return(u)
}

get_wts = function(Xs, mu){
  wts = f(Xs) / g(Xs, mu)
  return(wts)
}

get_ESS = function(weights){
  ess = (sum(weights)^2 / sum(weights^2))
  return(ess)
}



k_hat_of_theta <- function(theta, n) {
  Xs = sample_from_proposal(n, theta)
  wts = get_wts(Xs, theta)
  pareto_smooth(wts, M_keep = "default", return_k = TRUE)$k_hat
}


grad_k_hat_of_theta = function(theta, Xs, delta = sqrt(.Machine$double.eps)) {
  Xs_plus = Xs + delta
  
  wts = get_wts(Xs, theta)
  wts_plus = get_wts(Xs_plus, theta + delta)
  
  k_hat = pareto_smooth(wts, M_keep = "default", return_k = TRUE)$k_hat
  k_hat_plus = pareto_smooth(wts_plus, M_keep = "default", return_k = TRUE)$k_hat
  (k_hat_plus - k_hat) / delta
}

# # square root of machine epsilon
# delta = sqrt(.Machine$double.eps)
# 
# k_hat_of_theta(1, 1000)
# grad_k_hat_of_theta(1, 1000)
# 
# ran_grad_k_hat = function(theta, n, delta = 1e-2) {
#   k_hat = k_hat_of_theta(theta, n)
#   k_hat_plus = k_hat_of_theta(theta + delta, n)
#   
#   (k_hat_plus - k_hat) / delta
# }



update_proposal = function(x, mu, stepsize){
  grad = grad_k_hat_of_theta(mu, x)
  new_mu = mu - grad*stepsize
  
  return(new_mu)
}


set.seed(111)

#
# Initialize values
#

# Number of Monte Carlo samples
# MC = 100
MC = 1000

# Initial value for mu in proposal
mu = numeric(MC)
all_k_hats = numeric(MC)

mu[1] = 1.5

# let mu = 1

# Number of random sample at each iteration
# N = 100
N = 1000

xx = sample_from_proposal(N, mu[1])
wts = get_wts(xx, mu[1])
all_k_hats[1] = pareto_smooth(wts, M_keep = "default", return_k = TRUE)$k_hat

for(i in 2:MC){
  print(paste0(i, " out of ", MC))
  # You can change stepsize if you run into problem
  # but it should be ok with this one (William: I think so).
  # mu[i]  = update_proposal(xx, mu[i-1], 0.2*(1/i)^(0.55) )
  mu[i]  = update_proposal(xx, mu[i-1], 1*(1/i)^(0.90) )
  xx     = sample_from_proposal(N, mu[i])
  wts = get_wts(xx, mu[i])
  all_k_hats[i] = pareto_smooth(wts, M_keep = "default", return_k = TRUE)$k_hat
}
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

