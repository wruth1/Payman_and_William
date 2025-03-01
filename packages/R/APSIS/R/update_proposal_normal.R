#' Updating the Proposal, Generating New Values for Mean Parameter
#'
#' @param n The desired sample size to generate at each iteration. This argument
#'   is being passed from \code{\link{optimise_proposal_normal}} function.
#' @param mu The current value of generated mu that will be updated. This
#'   argument is being passed from \code{\link{optimise_proposal_normal}}
#'   function.
#' @param f The target function being passed from
#'   \code{\link{optimise_proposal_normal}}.
#' @param step.size The value of step size in optimization part when updating
#'   generated mean.
#'
#' @returns
#' @export
#'
update_proposal_normal = function(n, mu, f, step.size){

  # Generate a random sample from Normal(mu,1) for the current value of mu
  xsample       <- rnorm(n = n, mean = mu)

  # Evaluate target distribution at generated sample
  eval_target   <- f(xsample)

  # Evaluate proposal distribution at generated sample with current value of mu
  eval_proposal <- dnorm(x = xsample, mean = mu)

  # Compute weights
  weight        <- eval_target / eval_proposal

  # Calculate gradient value
  gt            <- estimate_gradient_normal(x = xsample, mu = mu, w = weight)

  # Compute new value for mu (update proposal)
  new_mu <- mu - (gt * step.size)

  # Return
  return(new_mu)

}
