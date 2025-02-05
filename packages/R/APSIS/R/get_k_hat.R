#' Estimate k for a given value of theta and a given dataset
#'
#' @param theta
#' @param X
#'
#' @returns
#'
#' @examples
get_k_hat <- function(theta, x) {

  # Compute K hat in formula (4) paper Zhang and Stephens
  each_term <- log(1 - theta * x)
  k_hat <- -mean(each_term)

  # Return the results
  return(k_hat)
}
