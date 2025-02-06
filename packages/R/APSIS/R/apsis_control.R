#' Auxiliary Function for Controlling Parameters in APSIS Algorithm
#'
#' @param epsilon A numeric value. the tolerance being used to check the
#'   absolute differences between two generated parameters from target
#'   distribution and stop the simulation. The default value is 1e-3.
#' @param max.iter A positive numeric value. The maximum number of Monte Carlo
#'   simulation that you want the algorithm to take. The default value is 1000.
#'   The algorithm stops when either desired tolerance achieved or max number of
#'   iteration reached, whichever comes first.
#' @param step.size A positive numeric value for step size value when updating
#'   proposal. The best practical suggestion is to set stepsize as \frac{1}{t}
#'   where t is the iteration number and this is the default behaviour when
#'   step.size = NULL.
#'
#' @return A list of three.
#' @export
#'
apsis_control = function(epsilon = 1e-4, max.iter = 1000, step.size = NULL){

  # Check if epsilon is a numeric positive value
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop('value of epsilon must be > 0')

  # Check if max.iter is a numeric positive value
  if (!is.numeric(max.iter) || max.iter <= 0)
    stop('maximum number of iterations must be > 0')

  # Create a list with control parameters
  res <- list(epsilon = epsilon, max.iter = max.iter, step.size = step.size)

  # Return the list
  return(res)

}
