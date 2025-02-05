#' Title
#'
#' @param x
#'
#' @returns
#'
#' @examples
estimate_theta = function(x){

  # Compute 0.25 quantile of x
  X_star      <- quantile(x, 0.25)

  # Get the maximum of x
  X_max       <- max(x)

  # Compute M according to Zhang and Stephens
  M_theta     <- round(20 + sqrt(length(x)))

  # Call get_thetas function to compute thetas according to (5) in Zhang and
  # Stephens paper
  all_thetas  <- get_thetas(X_star, X_max, M_theta)

  # Call get_all_weights function to compute weights of each theta for formula
  # (6) in Zhang and Stephens paper
  all_weights <- get_all_weights(all_thetas, X)

  # Finally compute MLE of theta according to forumla (6)
  theta_hat   <- sum(all_weights * all_thetas)

  # Return results
  return(theta_hat)

}
