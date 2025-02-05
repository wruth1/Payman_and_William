#' Compute intermediate estimates of theta for fitting the Pareto distribution
#'
#' @param X_star
#' @param X_max
#' @param M_theta
#'
#' @returns
#'
#' @examples
get_thetas = function(X_star, X_max, M_theta){

  # Compute thetas for formula (6) in Zhang and Stephens paper
  # Note that we have the formula in 3 parts
  part1 <- 1 / X_max
  part2 <- 1 - sqrt(M_theta / (1:M_theta - 0.5))
  part3 <- 3 * X_star

  # Calculate all_thetas
  all_thetas <-  part1 + (part2 / part3)

  # Return results
  return(all_thetas)
}
