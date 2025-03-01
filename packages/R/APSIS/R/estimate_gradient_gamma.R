#' Title
#'
#' @param x
#' @param mu
#' @param w
#'
#' @returns
#' @export
#'
#' @l
estimate_gradient_gamma = function(x, l, w){

  # The calculation is based on section 3.2 of the paper ”Convergence rates for
  # optimised adaptive importance samplers by Akyildiz and Miguez, 2022 We
  # estimate the gradient of "effective sample size" by an average of generated
  # sample and weights
  temp <- ( x - (1/l) ) * w

  # Compute sample mean
  gt   <- mean(temp)

  # Return
  return(gt)
}
