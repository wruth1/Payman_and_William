

library(latex2exp)
library(pbapply)
library(parallel)
library(magrittr)
library(evmix)  # For generalized Pareto distribution functions (e.g. dgpd)

source("src/Pareto_Smoothing.R")


plot_dir = "Presentations/MacEwan 2025/Figures/"


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

get_log_wts = function(Xs, sigma){
  log_f = dnorm(Xs, mean = 0, sd = 1, log=T)
  log_g = dnorm(Xs, mean = 0, sd = sigma, log=T)

  log_wts = log_f - log_g
  return(log_wts)
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
#                                Running Example                               #
# ---------------------------------------------------------------------------- #
#! Old setting for presentation was: set.seed(1), theta_bad = 0.6, theta_good = 2
# set.seed(1)
set.seed(11)

n_EG = 1000

theta_good = 2
theta_bad = 1/theta_good
# theta_bad = 0.6

Xs_good = sample_from_proposal(n_EG, theta_good)
Xs_bad = sample_from_proposal(n_EG, theta_bad)

wts_good = get_wts(Xs_good, sigma = theta_good)
wts_bad = get_wts(Xs_bad, sigma = theta_bad)


(E_hat_good = mean(wts_good * Xs_good^2))
(E_hat_bad = mean(wts_bad * Xs_bad^2))

mean(wts_good * (Xs_good - E_hat_good)^2)
mean(wts_bad * (Xs_bad - E_hat_bad)^2)

#* Histograms of weights
pdf(paste0(plot_dir, "Wt Hist.pdf"), width=10, height=7)
par(mfrow=c(1,2))
hist(wts_good, main = TeX(r'($G_1 = N(0, 2^2)$)'), xlab = 'Weights', freq=F, breaks = 50)
hist(wts_bad, main = TeX(r'($G_2 = N(0, 0.5^5)$)'), xlab = 'Weights', freq=F, breaks = 50)
dev.off()
par(mfrow=c(1,1))


#* Effective Sample Size

get_ESS <- function(wts){
  return(sum(wts)^2 / sum(wts^2))
}

(ESS_good = get_ESS(wts_good))
(ESS_bad = get_ESS(wts_bad))



# ------------------------------- Truncated IS ------------------------------- #

#* Compute thresholds
thresh_good = sqrt(n_EG) * mean(wts_good)
thresh_bad = sqrt(n_EG) * mean(wts_bad)


#* Display threshold(s); only one is visible
pdf(paste0(plot_dir, "Wt Hist - Thresh.pdf"), width=10, height=7)
par(mfrow=c(1,2))
hist(wts_good, main = TeX(r'($G_1 = N(0, 2^2)$)'), xlab = 'Weights', freq=F, breaks = 50)
hist(wts_bad, main = TeX(r'($G_2 = N(0, 0.5^5)$)'), xlab = 'Weights', freq=F, breaks = 50)
abline(v = thresh_bad, col = 'red', lwd = 2)
dev.off()
par(mfrow=c(1,1))


#* Apply truncation
wts_good_trunc = wts_good
wts_good_trunc[wts_good > thresh_good] = thresh_good
wts_bad_trunc = wts_bad
wts_bad_trunc[wts_bad > thresh_bad] = thresh_bad

pdf(paste0(plot_dir, "Wt Hist - Trunc.pdf"), width=10, height=7)
par(mfrow=c(1,2))
hist(wts_good_trunc, main = TeX(r'($G_1 = N(0, 2^2)$)'), xlab = 'Weights', freq=F, breaks = 50)
hist(wts_bad_trunc, main = TeX(r'($G_2 = N(0, 0.5^5)$)'), xlab = 'Weights', freq=F, breaks = 50)
dev.off()
par(mfrow=c(1,1))


#* ESS after truncating
(ESS_good_trunc = get_ESS(wts_good_trunc))
(ESS_bad_trunc = get_ESS(wts_bad_trunc))


# ----------------------------------- PSIS ----------------------------------- #
source("src/Pareto_Smoothing.R")
M_keep <- get_M_keep(wts_good)


#* Get thresholds
wts_good_decr = order(wts_good, decreasing = TRUE)
ind_keep_good = wts_good_decr[1:M_keep]
wts_good_keep = wts_good[ind_keep_good]
wts_good_thresh = min(wts_good_keep)

wts_bad_decr = order(wts_bad, decreasing = TRUE)
ind_keep_bad = wts_bad_decr[1:M_keep]
wts_bad_keep = wts_bad[ind_keep_bad]
wts_bad_thresh = min(wts_bad_keep)

#* Plot Thresholds
pdf(paste0(plot_dir, "Wt Hist - Pareto Thresh.pdf"), width=10, height=7)
par(mfrow=c(1,2))
hist(wts_good, main = TeX(r'($G_1 = N(0, 2^2)$)'), xlab = 'Weights', freq=F, breaks = 50)
abline(v = wts_good_thresh, col = 'red', lwd = 2)

hist(wts_bad, main = TeX(r'($G_2 = N(0, 0.5^5)$)'), xlab = 'Weights', freq=F, breaks = 50)
abline(v = wts_bad_thresh, col = 'red', lwd = 2)
dev.off()
par(mfrow=c(1,1))


#* Fit GPD smoother
GPD_pars_good = fit_GPD(wts_good_keep)
GPD_pars_bad = fit_GPD(wts_bad_keep)

# Note: Density functions are NA below the threshold value, u
tail_dens_good = function(w) dgpd(w, u = GPD_pars_good$mu_hat, sigmau = GPD_pars_good$sigma_hat, xi = GPD_pars_good$k_hat, phiu = M_keep / n_EG)
tail_dens_bad = function(w) dgpd(w, u = GPD_pars_bad$mu_hat, sigmau = GPD_pars_bad$sigma_hat, xi = GPD_pars_bad$k_hat, phiu = M_keep / n_EG)


#* Plot smoothed densities
pdf(paste0(plot_dir, "Wt Hist - Pareto Dens.pdf"), width=10, height=7)
par(mfrow=c(1,2))
hist(wts_good, main = TeX(r'($G_1 = N(0, 2^2)$)'), xlab = 'Weights', freq=F, breaks = 50)
abline(v = wts_good_thresh, col = 'red', lwd = 2)
curve(tail_dens_good, from = wts_good_thresh+0.001, to = max(wts_good), col = 3, add = T, lwd = 2)

hist(wts_bad, main = TeX(r'($G_2 = N(0, 0.5^5)$)'), xlab = 'Weights', freq=F, breaks = 50)
abline(v = wts_bad_thresh, col = 'red', lwd = 2)
curve(tail_dens_bad, from = wts_bad_thresh+0.01, to = max(wts_bad), col = 3, add = T, lwd = 2)
dev.off()
par(mfrow=c(1,1))



# Zoom-in on tail
tail_dens_good_plot <- function(w) tail_dens_good(w) * n_EG / M_keep
tail_dens_bad_plot <- function(w) tail_dens_bad(w) * n_EG / M_keep

pdf(paste0(plot_dir, "Wt Hist - Pareto Dens Zoom.pdf"), width=10, height=7)
par(mfrow=c(1,2))
hist(wts_good_keep, main = TeX(r'($G_1 = N(0, 2^2)$)'), xlab = 'Weights', freq=F, xlim = c(wts_good_thresh, max(wts_good_keep)))
abline(v = wts_good_thresh, col = 'red', lwd = 2)
curve(tail_dens_good_plot, from = wts_good_thresh+0.0001, to = max(wts_good_keep), col = 3, add = T, lwd = 2)

hist(wts_bad_keep, main = TeX(r'($G_2 = N(0, 0.5^5)$)'), xlab = 'Weights', freq=F, xlim = c(wts_bad_thresh, max(wts_bad_keep)), breaks = 50)
abline(v = wts_bad_thresh, col = 'red', lwd = 2)
curve(tail_dens_bad_plot, from = wts_bad_thresh+0.01, to = max(wts_bad_keep), col = 3, add = T, lwd = 2)
dev.off()
par(mfrow=c(1,1))



#* Compute and plot pareto smoothed weights
wts_good_smooth = pareto_smooth(wts_good, M_keep)
wts_bad_smooth = pareto_smooth(wts_bad, M_keep)

pdf(paste0(plot_dir, "Wt Hist - Pareto Smooth.pdf"), width=14, height=7)
par(mfrow=c(1,2))
hist(wts_good_smooth, main = TeX(r'($G_1 = N(0, 2^2)$)'), xlab = 'Weights', freq=F, breaks = 50)

hist(wts_bad_smooth, main = TeX(r'($G_2 = N(0, 0.5^5)$)'), xlab = 'Weights', freq=F, breaks = 50)
dev.off()
par(mfrow=c(1,1))

# Zoom-in on tails
smoothed_good_keep = wts_good_smooth[ind_keep_good]
smoothed_bad_keep = wts_bad_smooth[ind_keep_bad]

pdf(paste0(plot_dir, "Wt Hist - Pareto Smooth Zoom.pdf"), width=10, height=7)
par(mfrow=c(1,2))
hist(smoothed_good_keep, main = TeX(r'($G_1 = N(0, 2^2)$)'), xlab = 'Weights', freq=F, xlim = c(wts_good_thresh, max(smoothed_good_keep)))
abline(v = wts_good_thresh, col = 'red', lwd = 2)
curve(tail_dens_good_plot, from = wts_good_thresh+0.0001, to = max(smoothed_good_keep), col = 3, add = T, lwd = 2)


hist(smoothed_bad_keep, main = TeX(r'($G_2 = N(0, 0.5^5)$)'), xlab = 'Weights', freq=F, xlim = c(wts_bad_thresh, max(wts_bad_keep)), breaks = 50)
abline(v = wts_bad_thresh, col = 'red', lwd = 2)
curve(tail_dens_bad_plot, from = wts_bad_thresh+0.01, to = max(wts_bad_keep), col = 3, add = T, lwd = 2)
dev.off()
par(mfrow=c(1,1))


#* ESS after Pareto smoothing
(ESS_smooth_good = get_ESS(wts_good_smooth))
(ESS_smooth_bad = get_ESS(wts_bad_smooth))


#* Compile all diagnostics
diagnostics_good = c(ESS_good, ESS_good_trunc, ESS_smooth_good, GPD_pars_good$k_hat)
diagnostics_bad = c(ESS_bad, ESS_bad_trunc, ESS_smooth_bad, GPD_pars_bad$k_hat)
both_diagnostics = cbind(diagnostics_good, diagnostics_bad)
colnames(both_diagnostics) = c("Good", "Bad")
rownames(both_diagnostics) = c("ESS", "ESS_trunc", "ESS_smooth", "k_hat")
print(both_diagnostics)



# ---------------------- Stochastic Approximation - ESS ---------------------- #

#! We optimize over zeta := log(sigma) so the support is unbounded

d_wt_sq_d_sigma <- function(sigma, x){
  A = (sigma^2 - x^2)/(sigma^3)
  B = get_wts(x, sigma)^2
  return(A * B)
}

#** Gradient of ESS wrt zeta
d_wt_sq_d_zeta <- function(zeta, x){
  sigma_grad = exp(zeta)
  wt_sq_grad = d_wt_sq_d_sigma(exp(zeta), x)

  return(wt_sq_grad * sigma_grad)
}

update_zeta_using_rho <- function(zeta_old, step_size, Xs, wts){
  some_grads = d_wt_sq_d_zeta(zeta_old, Xs)
  grad_rho_hat = mean(some_grads * wts)

  zeta_new = zeta_old - grad_rho_hat * step_size
  return(zeta_new)
}

update_sigma_using_rho <- function(sigma_old, step_size, Xs, wts){
  some_grads = d_wt_sq_d_sigma(sigma_old, Xs)
  grad_rho_hat = mean(some_grads * wts)

  sigma_new = sigma_old - grad_rho_hat * step_size
  return(sigma_new)
}




#* Run SA to maximize ESS (i.e. minimize rho, the expected squared weight)

set.seed(1)

zeta_old = log(0.9)                         #! Check if I can get away with log(0.5) or log(0.3). rho is infinite for sigma^2 < 0.5

step_size_power = 1

B_EG_ESS = 100

all_zetas = numeric(B_EG_ESS + 1)
all_zetas[1] = zeta_old

for(i in 1:B_EG_ESS){
  if(i %% 10 == 0) print(paste0(i, " of ", B_EG_ESS))

  sigma_old = exp(zeta_old)
  Xs = sample_from_proposal(n_EG, sigma_old)
  log_wts = get_log_wts(Xs, sigma_old)
  wts = get_wts(Xs, sigma_old)
  step_size = 1/i^step_size_power
  zeta_new = update_zeta_using_rho(zeta_old, step_size, Xs, wts)

  all_zetas[i+1] = zeta_new
  zeta_old = zeta_new

}

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


