

pacman::p_load(tidyr)

#
# Define the target distribution, f(x) standard Normal distribution.
# The function Takes x and returns p.d.f of N(0,1)
#
f = function(x){
  dnorm(x)
}


#
# Density of the proposal distribution.
# T with df degrees of freedom.
g = function(x, df){
    dt(x, df)
}

sim_G = function(N, df){
    rt(N, df)
}

sim_G_from_U = function(U, df){
  qt(U, df)
}

get_L1_ESS = function(W){
    N = length(W)

    W_p = subset(W, W>1/N)

    S_p = sum(W_p)
    N_p = length(W_p)

    return(-N * S_p + N_p + N)
}

# ToDo: Implement remaining norms
ESS_from_W = function(some_Ws, flavour = "L2"){
  if(flavour == "L2"){
    return(1/sum(some_Ws^2))
  } else if(flavour == "LInfty"){
    return(1 / max(some_Ws))
  } else if(flavour == "L1"){
    return(get_L1_ESS(some_Ws))
  } else if(flavour == "entropy"){
    return(exp(-sum(some_Ws * log(some_Ws))))
  }
}


ESS_from_X = function(df, some_Xs, flavour = "L2"){
    some_Ws = f(some_Xs) / g(some_Xs, df)

    some_Ws = some_Ws / sum(some_Ws)        #!!!!!!!!!!! Apply self-normalization

    return(ESS_from_W(some_Ws, flavour))
}


ESS_from_U = function(df, some_Us, flavour = "L2"){
  some_Xs = sim_G_from_U(some_Us, df)
  ESS_from_X(df, some_Xs, flavour)
}


ESS_from_N = function(df, N, flavour = "L2"){
    some_Xs = sim_G(N, df)
    ESS_from_X(df, some_Xs, flavour)
}



# ToDo: Analyze whether the common random numbers approach is appropriate for this calculation
# c: step size for finite difference. Denominator is 2c
fin_diff_grad = function(df, N, c, flavour = "L2"){

    some_Us = runif(N)

    A = ESS_from_U(df + c, some_Us, flavour = flavour)
    B = ESS_from_U(df - c, some_Us, flavour = flavour)

    return((A - B)/(2*c))
}

# c = get_c(2)
# fin_diff_grad(2, N, get_c(2))

#
# Define a function to update the proposal at each iteration, let's call this function update_proposal.
# This function takes the sample (x) and parameter (df) at each iteration and returns
# a new value for the parameter (df).
#
# (Some notes on stepsize:
# the best practical suggestion is to set stepsize as \frac{1}{t} where t is the iteration number.
# Usually the range for stepsize is between \frac{1}{sqrt(t)} and \frac{1}{t}
# Anything in this ranges should work it is a matter of time for convergence.)
#
# a: step size multiplier on gradient
# c: step size for finite difference. Denominator is 2c
update_proposal = function(df, N, a, c, flavour = "L2"){
  grad_hat = fin_diff_grad(df, N, c, flavour = flavour)
  
  df_new = df + a * grad_hat
  return(df_new)
}





# # Kiefer-Wolfowitz
# get_a = function(k) return(k^(-1))
# get_c = function(k) return(0.001 * k^(-1/3))

# # Aggressive
# get_a = function(k) return(k^(-0.8))
# get_c = function(k) return(k^(-0.25))

# Chen et al. (2024)
get_a = function(k, alpha = -0.501, a0 = 1) return(a0* k^(alpha))
get_c = function(k, alpha = -0.501, c0 = 0.001) return(c0 * k^(alpha))




# Number of SA iterations
MC = 1000
# MC = 1000

# Number of random sample at each iteration
# N = 100
N = 1000

#? If df ever lands outside a cmpt interval, project to the endpoints
df_min = 1   #! Might be better to use 1.5 to avoid Cauchy degeneracies
df_max = 200 #! Might need a larger value, as df -> infinity gives a normal


# flavour = "L1"
flavour = "L2"
# flavour = "LInfty"
# flavour = "entropy"


set.seed(111)

# Initial value for df in proposal
df = numeric(MC)
df[1] = 20



#
# sigma_0: Initial value

run_SA = function(theta_0, flavour,
                    theta_min = -Inf, theta_max = Inf,
                    N = 1000, MC = 500, 
                    a0 = 1, c0 = 0.001, 
                    alpha_a = -0.501, alpha_c = -0.501,
                    verbose = FALSE){

  # Initial value for theta in proposal
  theta = numeric(MC)
  theta[1] = theta_0


  for(i in 2:MC){
    if(verbose){
      if((i %% 50) == 0)  print(paste0(i, " out of ", MC))
      # print(paste0(i, " out of ", MC))
    }
      
      # Run update
      # theta[i]  = update_proposal(theta[i-1], N, get_a(i, alpha = -1, a0 = 10), get_c(i), flavour = flavour)
      theta_new = update_proposal(theta[i-1], N, get_a(i, alpha = alpha_a, a0 = a0), get_c(i, alpha = alpha_c, c0 = c0), flavour = flavour)

        if(theta_new < theta_min){
          theta[i] = theta_min
      } else if(theta_new > theta_max){
          theta[i] = theta_max
      } else{
          theta[i] = theta_new
      }

      # Check bounds
      
  }

  return(theta)
}



# df_L1 = df
# df_L2 = df
# df_LInfty = df
# df_entropy = df


make_trajectories = function(theta_0 = 20, theta_min = 1, theta_max = 200, N=1000, MC = 1000, a0 = 1, c0 = 0.001, alpha_a = -0.501, alpha_c = -0.501){


  print("L1")
  theta_L1 = run_SA(theta_0 = theta_0, flavour = "L1",
                      theta_min = theta_min, theta_max = theta_max,
                      N = N, MC = MC, 
                      a0 = a0, c0 = c0, 
                      alpha_a = alpha_a, alpha_c = alpha_c)


  print("L2")
  theta_L2 = run_SA(theta_0 = theta_0, flavour = "L2",
                      theta_min = theta_min, theta_max = theta_max,
                      N = N, MC = MC, 
                      a0 = a0, c0 = c0, 
                      alpha_a = alpha_a, alpha_c = alpha_c)



  print("LInfty")
  theta_LInfty = run_SA(theta_0 = theta_0, flavour = "LInfty",
                      theta_min = theta_min, theta_max = theta_max,
                      N = N, MC = MC, 
                      a0 = a0, c0 = c0, 
                      alpha_a = alpha_a, alpha_c = alpha_c)



  print("entropy")
  theta_entropy = run_SA(theta_0 = theta_0, flavour = "entropy",
                      theta_min = theta_min, theta_max = theta_max,
                      N = N, MC = MC, 
                      a0 = a0, c0 = c0, 
                      alpha_a = alpha_a, alpha_c = alpha_c)



  data_theta = tibble(i = 1:MC, L1 = theta_L1, L2 = theta_L2, Infty = theta_LInfty, entropy = theta_entropy) %>%
      pivot_longer(2:5, names_to = "flavour", values_to = "theta")

    return(data_theta)
}

data_theta = make_trajectories(MC = 200, a0 = 100)

ggplot(data_theta, aes(x = i, y = theta)) + geom_line() + facet_wrap(~flavour)



df_ave_L1 = cumsum(df_L1) / seq_along(df_L1)
df_ave_L2 = cumsum(df_L2) / seq_along(df_L2)
df_ave_LInfty = cumsum(df_LInfty) / seq_along(df_LInfty)
df_ave_entropy = cumsum(df_entropy) / seq_along(df_entropy)


    pivot_longer(2:5, names_to = "flavour", values_to = "df")

# ToDo: Make a plot of the true ESS as a function of sigma. Compute the true value by Monte Carlo with high precision




sigma_ave_L1 = cumsum(sigma_L1) / seq_along(sigma_L1)
sigma_ave_L2 = cumsum(sigma_L2) / seq_along(sigma_L2)
sigma_ave_LInfty = cumsum(sigma_LInfty) / seq_along(sigma_LInfty)
sigma_ave_entropy = cumsum(sigma_entropy) / seq_along(sigma_entropy)

data_sigma_ave = tibble(i = 1:MC, L1 = sigma_ave_L1, L2 = sigma_ave_L2, Infty = sigma_ave_LInfty, entropy = sigma_ave_entropy) %>%
    pivot_longer(2:5, names_to = "flavour", values_to = "sigma_ave")

ggplot(data_sigma_ave, aes(x = i, y = sigma_ave)) + geom_line() + facet_wrap(~flavour) + geom_hline(yintercept = 1)



