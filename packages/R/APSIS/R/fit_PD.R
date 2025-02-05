#' Title
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
fit_PD = function(x){

  theta_hat <- estimate_theta(x)
  k_hat     <- get_k_hat(theta_hat, x)
  sigma_hat <- k_hat / theta_hat

  res <- list(sigma_hat = sigma_hat, k_hat = k_hat)
  return(res)

}
