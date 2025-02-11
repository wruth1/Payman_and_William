

library(latex2exp)
library(pbapply)
library(parallel)
library(magrittr)

source("src/Pareto_Smoothing.R")


plot_dir = "Presentations/RichCon 2024/Figures/"


# ---------------------------------------------------------------------------- #
#                             Densities and Weights                            #
# ---------------------------------------------------------------------------- #

#
# Define the target distribution, f(x) standard Normal distribution.
# The function Takes x and returns p.d.f of N(0,1)
#
f = function(x, mu=0){
  ( 1/sqrt(2*pi) ) * exp( -0.5*(x - mu)^2 )
}


#
# Define the proposal function, g(x,mu) N(\mu,1)
# A function of x (sample) and mu parameter (mean)
# For any values of x and \mu, returns the p.d.f of N(\mu,1)
#
g = function(x, sigma){
  (1/(sigma * sqrt(2*pi))) * exp(-0.5 * (x)^2 / sigma^2)
}

# sigma derivative of log(g)
g_score = function(x, sigma){
  A = x^2 / sigma^3
  B = 1 / sigma
  return(A - B)
}


sample_from_proposal = function(n, sigma){
  u = rnorm(n = n, mean = 0, sd = sigma)
  return(u)
}

get_wts = function(Xs, sigma){
  wts = f(Xs) / g(Xs, sigma)
  return(wts)
}





# ---------------------------------------------------------------------------- #
#                          Pareto Smoothing Functions                          #
# ---------------------------------------------------------------------------- #




theta_2_k_hat <- function(theta, n) {
  Xs = sample_from_proposal(n, theta)
  wts = get_wts(Xs, theta)
  pareto_smooth(wts, M_keep = "default", return_k = TRUE)$k_hat
}

Xs_2_k_hat <- function(theta, Xs) {
  wts = get_wts(Xs, theta)
  pareto_smooth(wts, M_keep = "default", return_k = TRUE)$k_hat
}


theta_2_grad_k_hat = function(theta, n, delta = sqrt(.Machine$double.eps)) {
 
  theta_plus = theta + delta
  theta_minus = theta - delta

  Xs = sample_from_proposal(n, theta)  #! This isn't right, but it's a good enough approximation for now
  Xs_plus = Xs
  Xs_minus = Xs

  wts_plus = get_wts(Xs_plus, theta_plus)
  wts_minus = get_wts(Xs_minus, theta_minus)
  
  k_hat_plus = pareto_smooth(wts_plus, M_keep = "default", return_k = TRUE)$k_hat
  k_hat_minus = pareto_smooth(wts_minus, M_keep = "default", return_k = TRUE)$k_hat
  grad_hat = (k_hat_plus - k_hat_minus) / (2 * delta)

  return(grad_hat)
}

Xs_2_grad_k_hat = function(theta, Xs, delta = sqrt(.Machine$double.eps)) {
 
  Xs_plus = Xs                                                                          #! Need to more explicitly handle dependance of Xs on theta
  Xs_minus = Xs

  theta_plus = theta + delta
  theta_minus = theta - delta

  wts_plus = get_wts(Xs_plus, theta_plus)
  wts_minus = get_wts(Xs_minus, theta_minus)
  
  k_hat_plus = pareto_smooth(wts_plus, M_keep = "default", return_k = TRUE)$k_hat
  k_hat_minus = pareto_smooth(wts_minus, M_keep = "default", return_k = TRUE)$k_hat
  grad_hat = (k_hat_plus - k_hat_minus) / (2 * delta)

  return(grad_hat)
}

# ------------------------------- Test Gradient ------------------------------ #

#* Note: I tried theta = 2, 0.5 and 3. Much better results with 0.5
set.seed(1)
theta = 2
n_test = 10000

Xs_test = sample_from_proposal(n_test, theta)

(math_grad = Xs_2_grad_k_hat(theta, Xs_test, delta = (.Machine$double.eps)^(1/4)))
(num_grad = numDeriv::grad(func = Xs_2_k_hat, x = theta, Xs = Xs_test))

rel_err = abs((num_grad - math_grad) / num_grad)
testthat::expect_true(rel_err < 1e-2)



#! ------------- Account for dependence of X on theta in gradient ------------- #

n_test = 1000
B_test = 100000

some_full_k_hat_grads = numeric(B_test)

set.seed(1)




# Setup cluster
# cl = makeCluster(detectCores() - 2)
# cl = makeCluster(15)
cl = makeCluster(10)
# clusterExport(cl, c("N", "b_Y", "theta_Y", "b_M", "theta_M", "which_REs"))
clusterExport(cl, c("n_test", "theta", "Xs_2_k_hat", "Xs_2_grad_k_hat", "sample_from_proposal", "g_score", "get_wts", "f", "g"))
clusterEvalQ(cl, {
    source("src/Pareto_Smoothing.R")
})
clusterSetRNGStream(cl = cl, 123)



some_full_k_hat_grads = pbsapply(1:B_test, function(i){
  new_Xs = sample_from_proposal(n_test, theta)

  k_hat_grad = Xs_2_grad_k_hat(theta, new_Xs, delta = (.Machine$double.eps)^(1/4))

  k_hat = Xs_2_k_hat(theta, new_Xs)
  L_prime = sum(g_score(new_Xs, theta))

  this_full_k_hat_grad = k_hat_grad + L_prime * k_hat
  return(this_full_k_hat_grad)
}, cl = cl)

stopCluster(cl)

mean(some_full_k_hat_grads)
sd(some_full_k_hat_grads) / sqrt(B_test)

mean(some_full_k_hat_grads[1:(B_test/10)])
mean(some_full_k_hat_grads[(1 + B_test/10):(2*B_test/10)])


set.seed(1)
this_B = 1000
this_theta = 0.1
this_n = 1000
some_metrics = pbsapply(1:this_B, function(i){
  some_Xs = sample_from_proposal(this_n, this_theta)
  some_wts = get_wts(some_Xs, sigma = this_theta)
  this_rho = mean(some_wts^2) 
  this_k_hat = Xs_2_k_hat(this_theta, some_Xs)

  return(c(rho_hat = this_rho, k_hat = this_k_hat))
})

some_metric_summaries = apply(some_metrics, 1, function(x){
  return(c(mean = mean(x), SE = sd(x) / sqrt(this_B), COV = sd(x) / mean(x)))
}) 

print(some_metric_summaries)


# mean(some_wts^2) / (mean(some_wts)^2)



# #! Focuses on using Xs generated with the correct theta. Future problem
# set.seed(1)
# theta = 3
# n_test = 10000
# delta = (.Machine$double.eps)^(1/4)


# Xs = sample_from_proposal(n_test, theta)
# Xs_plus = Xs * (theta + delta) / theta
# Xs_minus = Xs * (theta - delta) / theta

# math_grad = Xs_2_grad_k_hat(theta, Xs_plus, Xs_minus, delta = delta)
# num_grad = numDeriv::grad(func = Xs_2_k_hat, x = theta, Xs = Xs)

# rel_err = abs((num_grad - math_grad) / num_grad)
# testthat::expect_true(rel_err < 1e-2)


# --------------- Functions involving log(sigma), denoted zeta --------------- #

zeta_sample_from_proposal = function(n, zeta){
  sigma = exp(zeta)
  Xs = sample_from_proposal(n, sigma)
  return(Xs)
}


zeta_2_k_hat <- function(zeta, n) {
  sigma = exp(zeta)
  theta_2_k_hat(sigma, n)
}

zeta_Xs_2_k_hat <- function(zeta, Xs) {
  sigma = exp(zeta)
  Xs_2_k_hat(sigma, Xs)
}

zeta_2_grad_k_hat <- function(zeta, n, delta = sqrt(.Machine$double.eps)) {
  sigma = exp(zeta)
  theta_2_grad_k_hat(sigma, n, delta)
}

zeta_Xs_2_grad_k_hat <- function(zeta, Xs, delta = sqrt(.Machine$double.eps)) {
  sigma = exp(zeta)
  d_k_d_sigma = Xs_2_grad_k_hat(sigma, Xs, delta)

  d_sigma_d_zeta = exp(zeta)
  grad_hat = d_k_d_sigma * d_sigma_d_zeta

  return(grad_hat)
}

# ---------------------------- Test zeta gradient ---------------------------- #


set.seed(1)
zeta = 1.5
n_test = 10000

Xs = zeta_sample_from_proposal(n_test, zeta)

(math_grad = zeta_Xs_2_grad_k_hat(zeta, Xs, delta = (.Machine$double.eps)^(1/3)))
(num_grad = numDeriv::grad(func = zeta_Xs_2_k_hat, x = zeta, Xs = Xs))

(rel_err = abs((num_grad - math_grad) / num_grad))
testthat::expect_true(rel_err < 1e-2)



# ---------------------------------------------------------------------------- #
#                           Stochastic Approximation                           #
# ---------------------------------------------------------------------------- #

n_SA = 1000

step_size_power = 0.75   # Step size at iteration i is 1/i^step_size_power
#* Must be >= 0.75
grad_step_size_power = 0.25   # Radius of interval whose endpoints are used for finite difference approximation
#* Alt: 0.25

# Test that step size powers are valid
if(step_size_power > 1) stop("Step sizes must have divergent sum (i.e. step size power must be at most 1)")
if(step_size_power + grad_step_size_power <= 1) stop("Product of step sizes must have convergent sum (i.e. step size powers must satisfy step_size_power + grad_step_size_power > 1)")
if(step_size_power - grad_step_size_power <= 0.5) stop("Ratio of squared step sizes must have convergent sum (i.e. step size powers must satisfy step_size_power - grad_step_size_power > 0.5)")



# Setup algorithm
zeta_old = 3

num_iterations = 3000

zeta_traj = numeric(num_iterations+1)
zeta_traj[1] = zeta_old


# Xs = zeta_sample_from_proposal(n_SA, zeta_old)

k_hat_traj = numeric(num_iterations)
# k_hat_traj[1] = zeta_Xs_2_k_hat(zeta_old, Xs)

grad_k_hat_traj = numeric(num_iterations)
# grad_k_hat_traj[1] = zeta_Xs_2_grad_k_hat(zeta_old, Xs)


zeta_2_grad_k_hat(zeta_old, n_SA)


for(i in 1:num_iterations){
  if(i %% 100 == 0) print(paste0("Iteration ", i, " of ", num_iterations))
  
  this_Xs = zeta_sample_from_proposal(n_SA, zeta_old)
  this_grad_hat = zeta_Xs_2_grad_k_hat(zeta_old, this_Xs, delta = ((.Machine$double.eps)^(1/3))/(i^grad_step_size_power))


  # some_grad_hats_4 = replicate(100,   zeta_2_grad_k_hat(zeta_old, n_SA, delta = ((.Machine$double.eps)^(1/4))))
  # some_grad_hats_3 = replicate(100,   zeta_2_grad_k_hat(zeta_old, n_SA, delta = ((.Machine$double.eps)^(1/3))))
  # some_grad_hats_2 = replicate(100,   zeta_2_grad_k_hat(zeta_old, n_SA, delta = ((.Machine$double.eps)^(1/2))))

  # # zeta_2_grad_k_hat(zeta_old, n_SA, delta = ((.Machine$double.eps)^(1/4)))

  # mean(some_grad_hats_4)
  # mean(some_grad_hats_3)
  # mean(some_grad_hats_2)
  # sd(some_grad_hats_4) / 10
  # sd(some_grad_hats_3) / 10
  # sd(some_grad_hats_2) / 10
  
  this_step_len = 1/i^step_size_power
  
  zeta_new = zeta_old - this_step_len * this_grad_hat   #* Subtracting the gradient here moves us downhill



  k_hat_traj[i] = zeta_2_k_hat(zeta_new, n_SA)
  grad_k_hat_traj[i] = this_grad_hat

  zeta_traj[i+1] = zeta_new
  zeta_old = zeta_new
}

plot(zeta_traj)
# plot(zeta_traj[-1])
plot(exp(zeta_traj), ylab = "sigma")
# plot(exp(zeta_traj[-1]), ylab = "sigma")
plot(k_hat_traj)
plot(grad_k_hat_traj)

cummean <- function(X) cumsum(X) / (1:length(X))
plot(cummean(zeta_traj))


k_hat_traj



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

