

# This script computes the derivative of (1-F(w))/f(w) at a given value for w
# with both numerical and theoretical approach.
# The theoretical steps are shown in  Gamma LRT Appendix.pdf file in Notes folder.
# The numerical calculations are done by numDeriv package in R
# Some of the functions used in this script are shown in section 2.3 of Gamma LRT.pdf file.
# I put a comment in each function to refer to each pdf file.

# Load libraries for Lambert W function and numerical derivations
library(pracma)
library(numDeriv)

options(digits = 22)

# Set some of the parameters to constant values
lambda = 2
beta = 1
alpha = 3

gamma = 1/(1/beta - 1/lambda)
c = lambda / (beta^alpha * gamma(alpha))



# alpha <- 3
# c     <- 3
# gamma <- 3
A     <- (-1)/(gamma * (alpha-1) * c^(alpha-1))


w = 1
w_star = c * (gamma * (alpha - 1) * exp(-1))^(alpha - 1)


#? A note about the two branches of W:
#?      lambertWn <= lambertWp
#?      I expect that the one with a "p" is the principal branch. More importantly for us, this determines which is W_l and which is W_u. Just use alphabetical ordering.

#! Warning: lambertWp breaks for arguments too close to -1/e (it doesn't throw an error, it just runs forever). For us, that means we can't evaluate the survival function too close to the upper bound on W. I hope this won't be an issue. Otherwise, we may need to find a different implementation of W_u.

# Defined in formula 38 in section 2.2 of pdf file
Xlw = function(w, gamma, alpha, c){
  
  point <- ( -1/(gamma*(alpha-1)) ) * (w/c)^{1/(alpha-1)} 
  res   <- -gamma * (alpha - 1) * lambertWp(point)
  return(res)
  
}


# Defined in formula 38 in section 2.2 of pdf file
Xuw = function(w, gamma, alpha, c){
  
  point <- ( -1/(gamma*(alpha-1)) ) * (w/c)^{1/(alpha-1)} 
  res   <- -gamma * (alpha - 1) * lambertWn(point)
  return(res)
  
}

# Defined in formula 39 in section 2.3 of pdf file
SurvivalFunction = function(w, lambda, gamma, alpha, c){
  
  res <- exp( -(Xlw(w, gamma, alpha, c)/lambda) ) - exp( -(Xuw(w, gamma, alpha, c)/lambda) )
  return(res)
  
}



# Defined in formula 53 in section 2.3 of pdf file
f = function(w, lambda, A, gamma, alpha, c){
  
  eval_point <- A * w^{1 / (alpha-1)}
  term1 <- ( lambertWn(eval_point) / (1 + lambertWn(eval_point)) ) * exp( -(Xuw(w, gamma, alpha, c)/lambda) )
  term2 <- ( lambertWp(eval_point) / (1 + lambertWp(eval_point)) ) * exp( -(Xlw(w, gamma, alpha, c)/lambda) )
  output = gamma * (term1 - term2) / (lambda * w)
  return(output)
  
}

w = 1

some_Ws = seq(0.1, 2, by = 0.01)

num_grad = - grad(func = SurvivalFunction, x = some_Ws, lambda = lambda, gamma = gamma, alpha = alpha, c = c)
math_grad = f(some_Ws, lambda, A, gamma, alpha, c)

max((num_grad - math_grad) / num_grad)





# Defined in Gamma LRT Appendix.pdf
gl = function(w, A, gamma, alpha, c){
  
  eval_point <- A * w^{1 / (alpha-1)}
  num <- 1 - (alpha-1)*(1 + lambertWn(eval_point))  
  den <- (alpha - 1) * (1 + lambertWn(eval_point))^2
  return(num/den)
  
}

# Defined in Gamma LRT Appendix.pdf
gu = function(w, A, gamma, alpha, c){
  
  eval_point <- A * w^{1 / (alpha-1)}
  num <- 1 - (alpha-1)*(1 + lambertWp(eval_point))  
  den <- (alpha - 1) * (1 + lambertWp(eval_point))^2
  return(num/den)
  
}

# Defined in Gamma LRT Appendix.pdf
fprime = function(w, A, gamma, alpha, c){
  
  eval_point <- A * w^{1 / (alpha-1)}
  term1 <- (1/(w)^2) * (lambertWn(eval_point) / (1 + lambertWn(eval_point)) ) * exp( -(Xlw(w, gamma, alpha, c)/gamma) ) * gl(w, A, gamma, alpha, c)
  term2 <- (1/(w)^2) * (lambertWp(eval_point) / (1 + lambertWp(eval_point)) ) * exp( -(Xuw(w, gamma, alpha, c)/gamma) ) * gu(w, A, gamma, alpha, c)
  return(term1-term2)
}

# Defined in Gamma LRT Appendix.pdf
FromTheory = function(w, A, gamma, alpha, c){
  res = 1 + ( fprime(w, A, gamma, alpha, c)/f(w, A, gamma, alpha, c) ) * ( SurvivalFunction(w, gamma, alpha, c) / f(w, A, gamma, alpha, c) ) 
  return(-res)
}


FromNumerical = function(w, A, gamma, alpha, c){
  SurvivalFunction(w, gamma, alpha, c) / f(w, A, gamma, alpha, c)  
}


# The derivative of ( (1-F(w))/f(w) ) w.r.t w at point w = 1 according to the theoretical derivation is:
w = 1
FromTheory(w, A, gamma, alpha, c)

# Numerical check for derivative of ( (1-F(w))/f(w) ) w.r.t w at point w = 1 is:
grad(func = FromNumerical, x = w, A = A, gamma = gamma, alpha = alpha, c = c)


#
# The accuracy gets better for larger values of alpha
#
alpha <- 10
c     <- 3
gamma <- 3
A     <- (-1)/(gamma * (alpha-1) * c^{alpha-1})

FromTheory(w, A, gamma, alpha, c)
grad(func = FromNumerical, x = w, A = A, gamma = gamma, alpha = alpha, c = c)


#
# I am not happy with the results of some parameter combinations for example: 
#

alpha <- 2
c     <- 3
gamma <- 3
A     <- (-1)/(gamma * (alpha-1) * c^{alpha-1})

FromTheory(w, A, gamma, alpha, c)
grad(func = FromNumerical, x = w, A = A, gamma = gamma, alpha = alpha, c = c)
