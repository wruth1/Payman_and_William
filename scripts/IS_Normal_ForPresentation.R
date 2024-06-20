library(latex2exp)

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

set.seed(123)


pdf(paste0(plot_dir, "Wt Hist.pdf"), width=10, height=7)
par(mfrow=c(1,2))



# Proposal 1: N(0.1, 1)
x1 <- rnorm(n = M, mean = 0.1, sd = 1)
w1 <- p(x1) / q(x1, mu = 0.1)
#hist(w1, breaks = 10, main = "Weights for Proposal N(0.1,1)", xlab = "Weights")
hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)


# Proposal 2: N(2, 1)
x2 <- rnorm(n = M, mean = 2, sd = 1)
w2 <- p(x2) / q(x2, mu = 2)
hist(w2, breaks = 100, main = TeX(r'($G_2 = N(2,1)$)'), xlab = 'Weights', freq=F)

dev.off()

# Estimated expectations
E_hat1 = mean(x1^2 * w1)
E_hat2 = mean(x2^2 * w2)


ESS1 = sum(w1)^2 / sum(w1^2)
ESS2 = sum(w2)^2 / sum(w2^2)



# Truncated IS

thresh = sqrt(M)

## Plot histograms with threshold
pdf(paste0(plot_dir, "Wt Hist - Thresh.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)
hist(w2, breaks = 100, main = TeX(r'($G_2 = N(2,1)$)'), xlab = 'Weights', freq=F)
abline(v = thresh, col = 'red', lwd = 2)

dev.off()

## Plot histograms with truncated weights
## Use same x-axis scale as un-truncated weights
w2_trunc = w2
w2_trunc[w2_trunc > thresh] = thresh

pdf(paste0(plot_dir, "Wt Hist - Trunc.pdf"), width=10, height=7)
par(mfrow=c(1,2))

hist(w1, 'fd', main = TeX(r'($G_1 = N(0.1,1)$)'), xlab = 'Weights', freq=F)
hist(w2_trunc, breaks = 100, main = TeX(r'($G_2 = N(2,1)$)'), xlab = 'Weights', xlim = c(0, 60), freq=F)

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

hist(w2, 'fd', main = TeX(r'($G_2 = N(2,1)$)'), xlab = 'Weights', freq=F)
abline(v = w2_thresh, col = 'red', lwd = 2)

dev.off()



library(evmix)
