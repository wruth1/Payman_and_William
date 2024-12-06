# This script computes the derivative of (1-F(w))/f(w) at a given value for w
# with both numerical and theoretical approach.
# The theoretical steps are shown in  Gamma LRT Appendix.pdf file in Notes folder.
# The numerical calculations are done by numDeriv package in R
# Some of the functions used in this script are shown in section 2.3 of Gamma LRT.pdf file.
# I put a comment in each function to refer to each pdf file.

# Load libraries for Lambert W function and numerical derivations
library(pracma)
library(numDeriv)

options(digits = 18)

# Set some of the parameters to constant values
alpha <- 3
c     <- 3
gamma <- 3
A     <- (-1)/(gamma * (alpha-1) * c^{alpha-1})


# Defined in formula 38 in section 2.2 of pdf file
Xlw = function(w, gamma, alpha, c){
  
  point <- ( -1/(gamma*(alpha-1)) ) * (w/c)^{1/(alpha-1)} 
  res   <- gamma * (alpha - 1) * lambertWn(point)
  return(res)
  
}


# Defined in formula 38 in section 2.2 of pdf file
Xuw = function(w, gamma, alpha, c){
  
  point <- ( -1/(gamma*(alpha-1)) ) * (w/c)^{1/(alpha-1)} 
  res   <- gamma * (alpha - 1) * lambertWp(point)
  return(res)
  
}


# Defined in formula 39 in section 2.3 of pdf file
SurvivalFunction = function(w, gamma, alpha, c){
  
  res <- exp( -(Xuw(w, gamma, alpha, c)/gamma) ) - exp( -(Xlw(w, gamma, alpha, c)/gamma) )
  return(res)
  
}


# Defined in formula 53 in section 2.3 of pdf file
f = function(w, A, gamma, alpha, c){
  
  eval_point <- A * w^{1 / (alpha-1)}
  term1 <- ( lambertWn(eval_point) / (1 + lambertWn(eval_point)) ) * (1/w) * exp( -(Xlw(w, gamma, alpha, c)/gamma) )
  term2 <- ( lambertWp(eval_point) / (1 + lambertWp(eval_point)) ) * (1/w) * exp( -(Xuw(w, gamma, alpha, c)/gamma) )
  return(term1 - term2)
  
}

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


b = 2
lambda = 1
c = 



#
# Test pdf f(w)
#
w = 0
grad(func = SurvivalFunction, x = w, gamma = gamma, alpha = alpha, c = c, method.args=list(r=6) )

f(w,A,gamma,alpha,c)

SurvivalFunction(w=1e-20, gamma, alpha, c)
y <- SurvivalFunction(w=w, gamma, alpha, c)
plot(w,y)


gamma * c * (alpha-1)



