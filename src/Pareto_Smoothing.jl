
export fit_PD, fit_GPD
export get_all_GPD_quantiles
export pareto_smooth_all, pareto_smooth, pareto_smooth!

using Random
using Distributions
using Statistics
using LogExpFunctions




# Random.seed!(1)

# M = 100
# X_raw = rand(Uniform(0, 1), M)
# X = X_raw ./ sum(X_raw)




# X_star = quantile(X, 0.25)
# X_max = maximum(X)

# M_theta = Int(20 + floor(sqrt(M)))

# ---------------------------------------------------------------------------- #
#                      Fit Generalized Pareto Distribution                     #
# ---------------------------------------------------------------------------- #

"""
Compute intermediate estimates of theta for fitting the Pareto distribution.
"""
function get_thetas(X_star, X_max, M_theta)
    all_thetas = 1/X_max .+ [1 - sqrt(M_theta / (j - 0.5)) for j in 1:M_theta] ./ (3 * X_star)
    return all_thetas
end


"""
Estimate k for a given value of theta and a given dataset, X.
"""
function get_k_hat(theta, X)
    each_term = [log(1 - theta * this_X) for this_X in X]
    k_hat = - mean(each_term)
    return k_hat
end


"""
Compute the profile log-likelihood for theta.
"""
function pareto_prof_log_lik(theta, X)
    this_k_hat = get_k_hat(theta, X)
    
    A = log(theta / this_k_hat)
    B = this_k_hat - 1

    output = length(X) * (A + B)
    return output
end



"""
Normalize and exponentiate the given log weights.
"""
function normalize_weights(all_log_weights)
    log_norm_const = logsumexp(all_log_weights)
    return exp.(all_log_weights .- log_norm_const)
end


"""
Compute log-weight for each intermediate theta estimate.
"""
function get_all_log_weights(all_thetas, X)
    all_log_weights = [pareto_prof_log_lik(this_theta, X) for this_theta in all_thetas]
    return all_log_weights
end


"""
Compute weight for each intermediate theta estimate.
"""
function get_all_weights(all_thetas, X)
    all_log_weights = get_all_log_weights(all_thetas, X)
    all_weights = normalize_weights(all_log_weights)
    return all_weights
end

"""
Estimate theta for a given dataset, X, using the method of Zhang and Stephens.
"""
function estimate_theta(X)
    X_star = quantile(X, 0.25)
    X_max = maximum(X)
    M_theta = Int(20 + floor(sqrt(length(X))))
    all_thetas = get_thetas(X_star, X_max, M_theta)
    all_weights = get_all_weights(all_thetas, X)
    theta_hat = sum(all_weights .* all_thetas)
    return theta_hat
end


"""
Estimate parameters for the Pareto distribution using the method of Zhang and Stephens.
"""
function fit_PD(X)
    theta_hat = estimate_theta(X)
    k_hat = get_k_hat(theta_hat, X)
    sigma_hat = k_hat / theta_hat
    return sigma_hat, k_hat
end


"""
Estimate parameters for the generalized Pareto distribution of Vehtari et al.
Note: This consists of translating the data to start at zero, then fitting the PD of Zhang and Stephens, then multiplying their k_hat by -1.
"""
function fit_GPD(X)
    mu_hat = minimum(X)
    X_trans = X .- mu_hat
    sigma_hat, neg_k_hat = fit_PD(X_trans)
    k_hat = - neg_k_hat

    return sigma_hat, k_hat, mu_hat
end



# ---------------------------------------------------------------------------- #
#                           Compute Quantiles of GPD                           #
# ---------------------------------------------------------------------------- #

"""
Compute the j/M th quantile of the GPD with given parameters.
"""
function get_one_GPD_quantile(j, M, sigma, k, mu)
    A1 = (j - 0.5) / M
    A = (1 - A1)^(-k)

    B = sigma * (A - 1) / k

    output = mu + B
    return output
end

"""
Compute M equally spaced quantiles of the given GPD.
"""
function get_all_GPD_quantiles(M, sigma, k, mu, order="descend")
    all_quantiles = [get_one_GPD_quantile(j, M, sigma, k, mu) for j in 1:M]
    if order=="descend"
        reverse!(all_quantiles)
    end
    return all_quantiles
end



# ---------------------------------------------------------------------------- #
#                           Perform Pareto Smoothing                           #
# ---------------------------------------------------------------------------- #

"""
Regularize fitted k value by shrinking it toward 0.5.
"""
function PSIS_regularize(k_hat, M_keep)
    A = M_keep * k_hat
    B = 5
    C = 10 + M_keep

    return (A + B)/C
end

"""
Fit GPD to all of X and return its quantiles. Outputted quantiles are in descending order.
"""
function pareto_smooth_all(X, return_k=false, regularize=true)
    sigma_hat, k_hat, mu_hat = fit_GPD(X)

    # Check for unacceptably large k_hat
    if k_hat > 0.7
        @warn "k_hat = $k_hat is large. Consider using a different proposal distribution."
    end

    # Regularize k_hat by shrinking toward 0.5
    if regularize
        k_hat = PSIS_regularize(k_hat, length(X))
    end
    all_quantiles = get_all_GPD_quantiles(length(X), sigma_hat, k_hat, mu_hat)
    if return_k
        return all_quantiles, k_hat
    else
        return all_quantiles
    end
end




"""
Get the default number of weights to retain for pareto smoothing.
"""
function get_M_keep(X)
    M_sqrt = ceil(Int, 3*sqrt(length(X)))
    M_prop = Int(0.2 * length(X))
    
    M_keep = min(M_sqrt, M_prop)
    return M_keep
end



# ToDo: Create versions of `pareto_smooth` and `pareto_smooth!` that specify the argument type for M_keep (i.e. integer or string).

"""
Fit GPD to largest M_keep elements of X. Return a copy of X with largest elements replaced by their smoothed values.
M_keep can either be an integer or "default". In the latter case, it is set to the recommended value in Vehtari et al.
"""
function pareto_smooth(X; M_keep, return_k=false, regularize=true)
    
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


"""
Fit GPD to largest M_keep elements of X. Return a copy of X with largest elements replaced by their smoothed values.
M_keep can either be an integer or "default". In the latter case, it is set to the recommended value in Vehtari et al.
"""
function pareto_smooth!(X; M_keep, return_k=false, regularize=true)

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
    X[ind_keep] = X_keep_smooth


    if return_k
        return X, k_hat
    else
        return X
    end
end



