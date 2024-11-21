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
alpha <- 3
c     <- 3
gamma <- 3
A     <- (-1)/(gamma * (alpha-1) * c^{alpha-1})


# Defined in formula 38 in section 2.2 of pdf file
Xlw = function(w){
  
  point <- ( -1/(gamma*(alpha-1)) ) * (w/c)^{1/(alpha-1)} 
  res   <- gamma * (alpha - 1) * lambertWn(point)
  return(res)
  
}


# Defined in formula 38 in section 2.2 of pdf file
Xuw = function(w){
  
  point <- ( -1/(gamma*(alpha-1)) ) * (w/c)^{1/(alpha-1)} 
  res   <- gamma * (alpha - 1) * lambertWp(point)
  return(res)
  
}


# Defined in formula 39 in section 2.3 of pdf file
SurvivalFunction = function(w){
  
  res <- exp( -(Xuw(w)/gamma) ) - exp( -(Xlw(w)/gamma) )
  return(res)
  
}


# Defined in formula 53 in section 2.3 of pdf file
f = function(w){
  
  eval_point <- A * w^{1 / (alpha-1)}
  term1 <- ( lambertWn(eval_point) / (1 + lambertWn(eval_point)) ) * (1/w) * exp( -(Xlw(w)/gamma) )
  term2 <- ( lambertWp(eval_point) / (1 + lambertWp(eval_point)) ) * (1/w) * exp( -(Xuw(w)/gamma) )
  return(term1 - term2)
  
}

# Defined in Gamma LRT Appendix.pdf
gl = function(w){
  
  eval_point <- A * w^{1 / (alpha-1)}
  num <- 1 - (alpha-1)*(1 + lambertWn(eval_point))  
  den <- (alpha - 1) * (1 + lambertWn(eval_point))^2
  return(num/den)
  
}

# Defined in Gamma LRT Appendix.pdf
gu = function(w){
  
  eval_point <- A * w^{1 / (alpha-1)}
  num <- 1 - (alpha-1)*(1 + lambertWp(eval_point))  
  den <- (alpha - 1) * (1 + lambertWp(eval_point))^2
  return(num/den)
  
}

# Defined in Gamma LRT Appendix.pdf
fprime = function(w){
  
  eval_point <- A * w^{1 / (alpha-1)}
  term1 <- (1/(w)^2) * (lambertWn(eval_point) / (1 + lambertWn(eval_point)) ) * exp( -(Xlw(w)/gamma) ) * gl(w)
  term2 <- (1/(w)^2) * (lambertWp(eval_point) / (1 + lambertWp(eval_point)) ) * exp( -(Xuw(w)/gamma) ) * gu(w)
  return(term1-term2)
}

# Defined in Gamma LRT Appendix.pdf
FromTheory = function(w){
  res = 1 + ( fprime(w)/f(w) ) * ( SurvivalFunction(w) / f(w) ) 
  return(-res)
}


FromNumerical = function(w){
  SurvivalFunction(w) / f(w)  
}


# The derivative of ( (1-F(w))/f(w) ) w.r.t w at point w = 1 according to the theoretical derivation is:
w = 1
FromTheory(w)

# Numerical check for derivative of ( (1-F(w))/f(w) ) w.r.t w at point w = 1 is:
grad(func = FromNumerical, x = w)


#
# The accuracy gets better for larger values of alpha
#
alpha <- 10
c     <- 3
gamma <- 3
A     <- (-1)/(gamma * (alpha-1) * c^{alpha-1})

FromTheory(w)
grad(func = FromNumerical, x = w)


#
# I am not happy with the results of some parameter combinations for example: 
#

alpha <- 2
c     <- 3
gamma <- 3
A     <- (-1)/(gamma * (alpha-1) * c^{alpha-1})

FromTheory(w)
grad(func = FromNumerical, x = w)
