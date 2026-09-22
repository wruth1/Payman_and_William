

# ---------------------------------------------------------------------------- #
#                         Easy - Normal Location Family                        #
# ---------------------------------------------------------------------------- #

# set.seed(1)

source("src/ESS Versions.r")


sigma = 0.1


# Kiefer-Wolfowitz
get_a = function(k) return(k^(-1))
get_c = function(k) return(k^(-1/3))

# # Aggressive
# get_a = function(k) return(k^(-0.8))
# get_c = function(k) return(k^(-0.25))



# Stochastic approximation to optimize a function, f, which can only be evaluated as f+E, with E ~ N(0, sigma^2). 



make_E = function(sigma = 1){
    return(rnorm(1, 0, sigma))
}


f_det = function(theta){
    -theta^2
}

# Evaluate f with noise (i.e. f+E). E is either provided or, if NULL, generated on the spot.
# Optionally, also return E (with f, as a list).
f_ran = function(theta, E = NULL, sigma=1, return_E = FALSE){
    
    if(is.null(E)){
        E = make_E(sigma)
    }

    f_val = f_det(theta)

    f_obs = f_val + E

    if(!return_E){
        return(f_obs)
    } else{
        return(list(f = f_obs, E = E))
    }
    
}


# Apply a finite difference approximation to the gradient of f (with noise)
# Either supply noise terms or they are generated. In latter case, different errors are used for step up and step down evaluations.
fin_diff_f = function(theta, c, E_up = NULL, E_down = NULL, sigma = 1, return_E = FALSE){
    if(is.null(E_up)) E_up = make_E(sigma)
    if(is.null(E_down)) E_down = make_E(sigma)

    A = f_ran(theta + c, E_up)
    B = f_ran(theta - c, E_down)

    grad_hat = (A - B) / (2*c)

    if(!return_E){
        return(grad_hat)
    } else{
        return(list(grad_hat = grad_hat, E_up = E_up, E_down = E_down))
    }
}



# k: iteration number
# n: size of sample from proposal
update_theta = function(theta_old, k, sigma = 1){
    a_k = get_a(k)
    
    c_k = get_c(k)
    grad_hat = fin_diff_f(theta_old, c_k, sigma = sigma)

    theta_new = theta_old + a_k * grad_hat
    return(theta_new)
}


# set.seed(1)

K = 1000

theta_init = 1
theta_old = theta_init

theta_traj = rep(0, times = K)
theta_traj[1] = theta_old

for(k in 2:K){
    theta_next = update_theta(theta_old, k, sigma = sigma)

    theta_traj[k] = theta_next
    theta_old = theta_next
}

# theta_traj


# par(mfrow = c(1,2))

# plot(theta_traj[-1], type = "l")
# abline(h = 0)

# f_traj = f_det(theta_traj)
# plot(f_traj[-1])

# par(mfrow = c(1,1))




# ---------------------------------------------------------------------------- #
#                         Analytic Gradient with Noise                         #
# ---------------------------------------------------------------------------- #


# Match time-independent component of iterationwise SD from finite difference version
sigma_grad = sigma / sqrt(2)

# Robbins-Monro suggestion
get_a_ana = function(k) 1/k

# Match Kiefer-Wolfowitz termwise SD
# get_a_ana = function(k) 1/k^(2/3)






grad_f_det = function(theta) -2 * theta

grad_f_ran = function(theta, E = NULL, sigma=1, return_E = FALSE){
    
    if(is.null(E)){
        E = make_E(sigma)
    }

    grad_f_val = grad_f_det(theta)

    grad_f_obs = grad_f_val + E

    if(!return_E){
        return(grad_f_obs)
    } else{
        return(list(grad_f = grad_f_obs, E = E))
    }
    
}

# k: iteration number
# n: size of sample from proposal
update_theta_ana = function(theta_old, k, sigma = 1){
    a_k = get_a_ana(k)
    
    grad_hat = grad_f_ran(theta_old, sigma = sigma)

    theta_new = theta_old + a_k * grad_hat
    return(theta_new)
}


# set.seed(1)

K_ana = K

theta_old = theta_init

theta_traj_ana = rep(0, times = K_ana)
theta_traj_ana[1] = theta_old

for(k in 2:K_ana){
    theta_next = update_theta_ana(theta_old, k, sigma = sigma)

    theta_traj_ana[k] = theta_next
    theta_old = theta_next
}

# theta_traj

theta_min = -0.05
theta_max = 0.05

y_min = -0.002

y_min_log = -22
y_max_log = -5

par(mfrow = c(2,3))

plot(theta_traj[-1], type = "l", main = "theta (Kfr-Wlf)", ylim = c(theta_min, theta_max))
abline(h = 0)

f_traj = f_det(theta_traj)
f_traj_log = log(-f_det(theta_traj))
plot(f_traj[-1], type = "l", main = "f (Kfr-Wlf)", ylim = c(y_min, 0))
plot(f_traj_log[-1], type = "l", main = "log(-f) (Kfr-Wlf)", ylim = c(y_min_log, y_max_log))


plot(theta_traj_ana[-1], type = "l", main = "theta (Rob-Mon)", ylim = c(theta_min, theta_max))
abline(h = 0)

f_traj_ana = f_det(theta_traj_ana)
f_traj_ana_log = log(-f_det(theta_traj_ana))
plot(f_traj_ana[-1], type = "l", main = "f (Rob-Mon)", ylim = c(y_min, 0))
plot(f_traj_ana_log[-1], type = "l", main = "log(-f) (Rob-Mon)", ylim = c(y_min_log, y_max_log))

par(mfrow = c(1,1))

dev.off()


# ---------------------------------------------------------------------------- #
#                          Apply Cumulative Averaging                          #
# ---------------------------------------------------------------------------- #


theta_traj_ave = cumsum(theta_traj) / seq_along(theta_traj)
f_traj_ave = f_det(theta_traj_ave)
f_traj_ave_log = log(- f_traj_ave)

theta_traj_ana_ave = cumsum(theta_traj_ana) / seq_along(theta_traj_ana)
f_traj_ana_ave = f_det(theta_traj_ana_ave)
f_traj_ana_ave_log = log(- f_traj_ana_ave)




par(mfrow = c(2,3))

plot(theta_traj_ave[-1], type = "l", main = "theta (Kfr-Wlf) - Averaged", ylim = c(theta_min, theta_max))
abline(h = 0)

plot(f_traj_ave[-1], type = "l", main = "f (Kfr-Wlf) - Averaged", ylim = c(y_min, 0))
plot(f_traj_ave_log[-1], type = "l", main = "log(-f) (Kfr-Wlf) - Averaged", ylim = c(y_min_log, y_max_log))


plot(theta_traj_ana_ave[-1], type = "l", main = "theta (Rob-Mon) - Averaged", ylim = c(theta_min, theta_max))
abline(h = 0)

plot(f_traj_ana_ave[-1], type = "l", main = "f (Rob-Mon) - Averaged", ylim = c(y_min, 0))
plot(f_traj_ana_ave_log[-1], type = "l", main = "log(-f) (Rob-Mon) - Averaged", ylim = c(y_min_log, y_max_log))

par(mfrow = c(1,1))

