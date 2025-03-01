#' Updating the Proposal, Generating New Values for Scale Parameter
#'
#' @param n The desired sample size to generate at each iteration. This argument
#'   is being passed from \code{\link{optimise_proposal_gamma}} function.
#' @param mu The current value of generated lambda that will be updated. This
#'   argument is being passed from \code{\link{optimise_proposal_gamma}}
#'   function.
#' @param f The target function being passed from
#'   \code{\link{optimise_proposal_gamma}}.
#' @param step.size The value of step size in optimization part when updating
#'   generated lambda.
#'
#' @returns
#' @export
#'
update_proposal_gamma = function(n, lambda, f, a, step.size){

  # Generate a random sample from Gamma(alpha,lambda) for the current value of
  # lambda
  xsample       <- rgamma(n = n, shape = a, scale = lambda)

  # Evaluate target distribution at generated sample
  eval_target   <- f(xsample)

  # Evaluate proposal distribution at generated sample with current value of mu
  eval_proposal <- dgamma(x = xsample, shape = a, scale = lambda)

  # Compute weights
  weight        <- eval_target / eval_proposal
  weight        <- weight * eval_target

  # Calculate gradient value
  gt            <- estimate_gradient_gamma(x = xsample, l = lambda, w = weight)

  # Compute new value for mu (update proposal)
  new_lambda <- lambda - (gt * step.size)

  # Return
  return(new_lambda)

}
