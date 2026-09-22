

# Kiefer-Wolfowitz
get_a = function(k) return(k^(-1))
get_c = function(k) return(k^(-1/3))


phi = function(x) x

f_dens = function(x) dnorm(x)
g_dens = function(x, mu, sigma) dnorm(x, mu, sigma)

get_W = function(x, mu, sigma) f_dens(x) / g_dens(x, mu, sigma)


sim_G = function(n, mu, sigma) rnorm(n, mu, sigma)


# Evaluate f with noise (i.e. f+E). E is either provided or, if NULL, generated on the spot.
# Optionally, also return E (with f, as a list).
f_ran = function(){
    
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