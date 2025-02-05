#' Pareto Smoothed Importance Sampling with a Family of Normal Distribution For
#' Proposal
#'
#' @param target
#' @param n
#' @param control
#' @param mu
#'
#' @returns
#' @export
#'
#' @description
#'
#' @examples
optimise_proposal_normal = function(target, n, mu = c(-10,10), control = NULL){

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
  # Set the values for current.iter which is a counter in the while statement
  # and check.cond to compare with epsilon
  #
  current.iter <- 0
  check.cond   <- FALSE

  while( (current.iter < max.iter) | check.cond ){

    ###########################################################################
    # Step1. Generate a random sample from Normal(mu,1) for a grid values of mu
    ###########################################################################

    #### First create a gird of values for some values of mu. This has a length
    #### of 201
    temp <- seq(from = mu[1], to = mu[2], by = 0.1)
    len  <- length(temp)

    #### Second make a matrix with one column and 201 rows where each row is one
    #### value of mu
    mu.grid <- matrix(data = temp, nrow = len, ncol = 1)

    #### Third take a sample of size n from proposal and make a matrix with n
    #### rows and 201 columns
    sample.from.proposal <- apply(X = mu.grid, MARGIN = 1, FUN = rnorm, n = n)
    ###########################################################################



    ###########################################################################
    #Step2. Compute weights. For each mu, we have one sample, and for each
    #sample we have a vector of weights
    ###########################################################################
    weights <- matrix(NA, nrow = n, ncol = len)
    for(i in 1:len){
      eval_target   <- target(x = sample.from.proposal[,i], mu.grid[i])
      eval_proposal <- dnorm(x = sample.from.proposal[,i], mean = mu.grid[i], sd = 1)
      weights[,i]   <- eval_target / eval_proposal
    }
    ###########################################################################



    # Step3. Implement Pareto Smoothing method to replace extreme weights

    ###########################################################################

  }


}

