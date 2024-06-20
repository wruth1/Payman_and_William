
library(matrixStats)

#### Fit Pareto Distribution ####

# Compute intermediate estimates of theta for fitting the Pareto distribution
get_thetas <- function(X_star, X_max, M_theta) {
  all_thetas <- 1 / X_max + (1 - sqrt(M_theta / (1:M_theta - 0.5))) / (3 * X_star)
  return(all_thetas)
}

# Estimate k for a given value of theta and a given dataset, X
get_k_hat <- function(theta, X) {
  each_term <- log(1 - theta * X)
  k_hat <- -mean(each_term)
  return(k_hat)
}

# Compute the profile log-likelihood for theta
pareto_prof_log_lik <- function(theta, X) {
  this_k_hat <- get_k_hat(theta, X)
  A <- log(theta / this_k_hat)
  B <- this_k_hat - 1
  output <- length(X) * (A + B)
  return(output)
}

# Normalize and exponentiate the given log weights
normalize_weights <- function(all_log_weights) {
  log_norm_const <- logSumExp(all_log_weights)
  return(exp(all_log_weights - log_norm_const))
}

# Compute log-weight for each intermediate theta estimate
get_all_log_weights <- function(all_thetas, X) {
  all_log_weights <- sapply(all_thetas, function(this_theta) pareto_prof_log_lik(this_theta, X))
  return(all_log_weights)
}

# Compute weight for each intermediate theta estimate
get_all_weights <- function(all_thetas, X) {
  all_log_weights <- get_all_log_weights(all_thetas, X)
  all_weights <- normalize_weights(all_log_weights)
  return(all_weights)
}

# Estimate theta for a given dataset, X, using the method of Zhang and Stephens
estimate_theta <- function(X) {
  X_star <- quantile(X, 0.25)
  X_max <- max(X)
  M_theta <- round(20 + sqrt(length(X)))
  all_thetas <- get_thetas(X_star, X_max, M_theta)
  all_weights <- get_all_weights(all_thetas, X)
  theta_hat <- sum(all_weights * all_thetas)
  return(theta_hat)
}

# Estimate parameters for the Pareto distribution using the method of Zhang and Stephens
fit_PD <- function(X) {
  theta_hat <- estimate_theta(X)
  k_hat <- get_k_hat(theta_hat, X)
  sigma_hat <- k_hat / theta_hat
  return(list(sigma_hat = sigma_hat, k_hat = k_hat))
}

# Estimate parameters for the generalized Pareto distribution of Vehtari et al.
fit_GPD <- function(X) {
  mu_hat <- min(X)
  X_trans <- X - mu_hat
  result <- fit_PD(X_trans)
  sigma_hat <- result$sigma_hat
  neg_k_hat <- result$k_hat
  k_hat <- -neg_k_hat
  return(list(sigma_hat = sigma_hat, k_hat = k_hat, mu_hat = mu_hat))
}



#### Compute Quantiles of GPD ####

# Compute the j/M th quantile of the GPD with given parameters
get_one_GPD_quantile <- function(j, M, sigma, k, mu) {
  A1 <- (j - 0.5) / M
  A <- (1 - A1)^(-k)
  B <- sigma * (A - 1) / k
  output <- mu + B
  return(output)
}

# Compute M equally spaced quantiles of the given GPD
get_all_GPD_quantiles <- function(M, sigma, k, mu, order = "descend") {
  all_quantiles <- sapply(1:M, function(j) get_one_GPD_quantile(j, M, sigma, k, mu))
  if (order == "descend") {
    all_quantiles <- rev(all_quantiles)
  }
  return(all_quantiles)
}

# Regularize fitted k value by shrinking it toward 0.5
PSIS_regularize <- function(k_hat, M_keep) {
  A <- M_keep * k_hat
  B <- 5
  C <- 10 + M_keep
  return((A + B) / C)
}

# Fit GPD to all of X and return its quantiles. Outputted quantiles are in descending order
pareto_smooth_all <- function(X, return_k = FALSE, regularize = TRUE) {
  fit <- fit_GPD(X)
  sigma_hat <- fit$sigma_hat
  k_hat <- fit$k_hat
  mu_hat <- fit$mu_hat
  
  # Check for unacceptably large k_hat
  if (k_hat > 0.7) {
    warning(paste0("k_hat = ", k_hat, " is large. Consider using a different proposal distribution."))
  }
  
  # Regularize k_hat by shrinking toward 0.5
  if (regularize) {
    k_hat <- PSIS_regularize(k_hat, length(X))
  }
  all_quantiles <- get_all_GPD_quantiles(length(X), sigma_hat, k_hat, mu_hat)
  if (return_k) {
    return(list(all_quantiles = all_quantiles, k_hat = k_hat))
  } else {
    return(all_quantiles)
  }
}

# Get the default number of weights to retain for pareto smoothing
get_M_keep <- function(X) {
  M_sqrt <- ceiling(3 * sqrt(length(X)))
  M_prop <- floor(0.2 * length(X))
  M_keep <- min(M_sqrt, M_prop)
  return(M_keep)
}

# Fit GPD to largest M_keep elements of X. Return a copy of X with largest elements replaced by their smoothed values
pareto_smooth <- function(X, M_keep = "default", return_k = FALSE, regularize = TRUE) {
  if (M_keep == "default") {
    M_keep <- get_M_keep(X)
  }
  
  # Extract largest elements of X for smoothing
  X_decr <- order(X, decreasing = TRUE)
  ind_keep <- X_decr[1:M_keep]
  X_keep <- X[ind_keep]
  
  # Apply smoothing
  result <- pareto_smooth_all(X_keep, TRUE, regularize)
  X_keep_smooth <- result$all_quantiles
  k_hat <- result$k_hat
  
  # Create smoothed copy of X
  X_smooth <- X
  X_smooth[ind_keep] <- X_keep_smooth
  
  for(i in seq_along(X_smooth)){
    if(X_smooth[i] > max(X)) X_smooth[i] <- max(X)
  }
  
  if (return_k) {
    return(list(X_smooth = X_smooth, k_hat = k_hat))
  } else {
    return(X_smooth)
  }
}

