#' Optimised Adaptive Importance Sampler with a Family of Normal Distribution for
#' Proposal Distribution
#'
#' @param target A target distribution. This must be a function of the sample
#'   (x) only. The `optimise_proposal_normal` function will evaluate this target
#'   function at randomly generated values of x. The target must return a
#'   numeric vector where each element is the evaluated value of x.
#' @param n The desired sample size to generate at each iteration.
#' @param mu0 The starting value for the mean parameter in the proposal normal
#'   distribution.
#' @param control A set of parameters to control algorithm convergence including
#'   stopping criteria and number of iterations.
#'
#' @returns
#' @export
#'
#' @description
#'
#' @examples
optimise_proposal_normal = function(target, n, mu0, control = NULL){

  #
  # Check if target is a numeric positive value
  #
  if( !is.function(target) ){
    stop('target must be a function')
  }

  #
  # Check if n is a numeric positive value
  #
  if (!is.numeric(n) || n < 2)
    stop('maximum number of iterations must be > 2')

  #
  # Set control parameters to default values
  #
  if( is.null(control) ){
    control <- apsis_control()
  }


  #
  # Initialize a vector to record simulated parameters of target distribution
  #
  param <- vector()


  #
  # Extract max.iter and epsilon from control list. 1) max.iter is the maximum
  # number of iteration that algorithm takes. In other words, the number of
  # iterations that we keep updating porposal. 2) epsilon is the tolerance that
  # we use to check if the difference between two consecutive generated
  # parameters is small enough
  #
  max.iter <- control$max.iter
  epsilon  <- control$epsilon


  #
  # Set the values before while loop:
  # t: keep track of iteration number
  # check: to evaluate convergence based on two consecutive mean values by
  # comparing to epsilon
  # Set the first element of param to initial value of mean set by user
  #
  t  <- 1
  check <- TRUE
  param[t] <- mu0


  while( (t < max.iter) | check ){

    # Update the proposal by finding the new value for mu
    param[t+1] <- update_proposal_normal(n = n, mu = param[t], f = target, step.size = 1/t)

    # Check if the convergence reached
    if( abs(param[t +1] - param[t]) < epsilon ){
      check <- FALSE
    }

    # Increase iteration
    t <- t + 1

  }

  # Return generated mu
  return(param)

}

