library(latex2exp)

set.seed(1)


plot_dir = "Presentations/RichCon 2024/Figures/"

M <- 1000

# Target distribution: N(0, 1)
p = function(x){
  return(dnorm(x))
}

# Proposal distribution: N(mu, 1)
q = function(x, mu){
  return(dnorm(x, mean = mu, sd = 1))
}



pdf(paste0(plot_dir, "Wt Hist.pdf"), width=10, height=7)
par(mfrow=c(1,2))



# Proposal 1: N(0.1, 1)
x1 <- rnorm(n = M, mean = 0.1, sd = 1)
w1 <- p(x1) / q(x1, mu = 0.1)
#hist(w1, breaks = 10, main = "Weights for Proposal N(0.1,1)", xlab = "Weights")
hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)


# Proposal 2: N(1.5, 1)
x2 <- rnorm(n = M, mean = 1.5, sd = 1)
w2 <- p(x2) / q(x2, mu = 1.5)
hist(w2, breaks = 100, main = TeX(r'($G_2 = N(1.5,1)$)'), xlab = 'Weights', freq=F)

dev.off()

# Estimated expectations
E_hat1 = mean(x1^2 * w1)
E_hat2 = mean(x2^2 * w2)


ESS1 = sum(w1)^2 / sum(w1^2)
ESS2 = sum(w2)^2 / sum(w2^2)



# Truncated IS

thresh1 = sqrt(M) * mean(w1)
thresh2 = sqrt(M) * mean(w2)

## Plot histograms with threshold
pdf(paste0(plot_dir, "Wt Hist - Thresh.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)
hist(w2, breaks = 100, main = TeX(r'($G_2 = N(1.5,1)$)'), xlab = 'Weights', freq=F)
abline(v = thresh2, col = 'red', lwd = 2)

dev.off()

## Plot histograms with truncated weights
## Use same x-axis scale as un-truncated weights
w2_trunc = w2
w2_trunc[w2_trunc > thresh2] = thresh2

pdf(paste0(plot_dir, "Wt Hist - Trunc.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)
hist(w2_trunc, breaks = 100, main = TeX(r'($G_2 = N(1.5,1)$)'), xlab = 'Weights', xlim = c(0, max(w2)), freq=F)

dev.off()


ESS2_trunc = sum(w2_trunc)^2 / sum(w2_trunc^2)





#### Pareto Smoothing ####

source("src/Pareto_Smoothing.R")
M_keep <- get_M_keep(w1)

w1_decr <- order(w1, decreasing = TRUE)
ind_keep1 <- w1_decr[1:M_keep]
w1_keep <- w1[ind_keep1]
w1_thresh = min(w1_keep)


w2_decr = order(w2, decreasing = TRUE)
ind_keep2 = w2_decr[1:M_keep]
w2_keep = w2[ind_keep2]
w2_thresh = min(w2_keep)



pdf(paste0(plot_dir, "Wt Hist - Pareto Thresh.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)
abline(v = w1_thresh, col = 'red', lwd = 2)

hist(w2, 'fd', main = TeX(r'($G_2 = N(1.5,1)$)'), xlab = 'Weights', freq=F)
abline(v = w2_thresh, col = 'red', lwd = 2)

dev.off()



library(evmix)

GPD_pars1 = fit_GPD(w1_keep)
GPD_pars2 = fit_GPD(w2_keep)

dens1 = function(w) dgpd(w, u = GPD_pars1$mu_hat, sigmau = GPD_pars1$sigma_hat, xi = GPD_pars1$k_hat, phiu = M_keep/M)
dens2 = function(w) dgpd(w, u = GPD_pars2$mu_hat, sigmau = GPD_pars2$sigma_hat, xi = GPD_pars2$k_hat, phiu = M_keep/M)


pdf(paste0(plot_dir, "Wt Hist - Pareto Dens.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)
abline(v = w1_thresh, col = 'red', lwd = 2)
curve(dens1, from = w1_thresh, to = max(w1), col = 3, add = T, lwd = 2)

hist(w2, 'fd', main = TeX(r'($G_2 = N(1.5,1)$)'), xlab = 'Weights', freq=F)
abline(v = w2_thresh, col = 'red', lwd = 2)
curve(dens2, from = w2_thresh, to = max(w2), col = 3, add = T, lwd = 2)

dev.off()


## Zoom-in on tails
pdf(paste0(plot_dir, "Wt Hist - Pareto Dens Zoom.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F, xlim = c(1.13, max(w1)), ylim = c(0,2))
abline(v = w1_thresh, col = 'red', lwd = 2)
curve(dens1, from = w1_thresh, to = max(w1), col = 3, add = T, lwd = 2)

hist(w2, 'fd', main = TeX(r'($G_2 = N(1.5,1)$)'), xlab = 'Weights', freq=F, xlim = c(4, max(w2)), ylim = c(0, 0.06))
abline(v = w2_thresh, col = 'red', lwd = 2)
curve(dens2, from = w2_thresh, to = max(w2), col = 3, add = T, lwd = 2)

dev.off()



## Pareto smoothed weights
w1_smooth = pareto_smooth(w1, M_keep)
w2_smooth = pareto_smooth(w2, M_keep)


pdf(paste0(plot_dir, "Wt Hist - Pareto Smooth.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1_smooth, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)
hist(w2_smooth, 'fd', main = TeX(r'($G_2 = N(1.5,1)$)'), xlab = 'Weights', freq=F)

dev.off()



### Zoom-in on tails
pdf(paste0(plot_dir, "Wt Hist - Pareto Smooth Zoom.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1_smooth, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F, xlim = c(1.13, max(w1_smooth)), ylim = c(0,2))

hist(w2_smooth, 'fd', main = TeX(r'($G_2 = N(1.5,1)$)'), xlab = 'Weights', freq=F, xlim = c(4, max(w2_smooth)), ylim = c(0, 0.06))
     
dev.off()


ESS1_smooth = sum(w1_smooth)^2 / sum(w1_smooth^2)
ESS2_smooth = sum(w2_smooth)^2 / sum(w2_smooth^2)
