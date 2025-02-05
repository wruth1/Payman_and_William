#' Apply Pareto Smoothing to Weights
#'
#' @param X
#' @param M_keep
#' @param return_k
#' @param regularize
#'
#' @returns
#' @export
#'
#' @examples
pareto_smoothing = function(X, M_keep, return_k = FALSE, regularize = TRUE){

  if M_keep == "default"
  M_keep = get_M_keep(X)
  end

  # Extract largest elements of X for smoothing
  X_decr = sortperm(X, rev=true)
  ind_keep = X_decr[1:M_keep]
  X_keep = X[ind_keep]

  # Apply smoothing
  X_keep_smooth, k_hat = pareto_smooth_all(X_keep, true, regularize)

  # Create smoothed copy of X
  X_smooth = deepcopy(X)
  X_smooth[ind_keep] = X_keep_smooth


  if return_k
  return X_smooth, k_hat
  else
    return X_smooth
  end
  end

}



