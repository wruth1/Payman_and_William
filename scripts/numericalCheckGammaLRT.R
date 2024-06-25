# Load library for Lambert W function and numerical derivations
library(pracma)
library(numDeriv)

options(digits = 22)

# Set some of the parameters to constant values
alpha <- 3
c     <- 3
gamma <- 3
A     <- (-1)/(gamma * (alpha-1) * c^{alpha-1})


Xlw = function(w){
  
  point <- ( -1/(gamma*(alpha-1)) ) * (w/c)^{1/(alpha-1)} 
  res   <- gamma * (alpha - 1) * lambertWn(point)
  return(res)
  
}


Xuw = function(w){
  
  point <- ( -1/(gamma*(alpha-1)) ) * (w/c)^{1/(alpha-1)} 
  res   <- gamma * (alpha - 1) * lambertWp(point)
  return(res)
  
}


OneMinusCapitalF = function(w){
  
  res <- exp( -(Xuw(w)/gamma) ) - exp( -(Xlw(w)/gamma) )
  return(res)
  
}

f = function(w){
  
  eval_point <- A * w^{1 / (alpha-1)}
  term1 <- ( lambertWn(eval_point) / (1 + lambertWn(eval_point)) ) * (1/w) * exp( -(Xlw(w)/gamma) )
  term2 <- ( lambertWp(eval_point) / (1 + lambertWp(eval_point)) ) * (1/w) * exp( -(Xuw(w)/gamma) )
  return(term1 - term2)
  
}

gl = function(w){
  
  eval_point <- A * w^{1 / (alpha-1)}
  num <- 1 - (alpha-1)*(1 + lambertWn(eval_point))  
  den <- (alpha - 1) * (1 + lambertWn(eval_point))^2
  return(num/den)
  
}


gu = function(w){
  
  eval_point <- A * w^{1 / (alpha-1)}
  num <- 1 - (alpha-1)*(1 + lambertWp(eval_point))  
  den <- (alpha - 1) * (1 + lambertWp(eval_point))^2
  return(num/den)
  
}

fprime = function(w){
  
  eval_point <- A * w^{1 / (alpha-1)}
  term1 <- (1/(w)^2) * (lambertWn(eval_point) / (1 + lambertWn(eval_point)) ) * exp( -(Xlw(w)/gamma) ) * gl(w)
  term2 <- (1/(w)^2) * (lambertWp(eval_point) / (1 + lambertWp(eval_point)) ) * exp( -(Xuw(w)/gamma) ) * gu(w)
  return(term1-term2)
}

FromTheory = function(w){
  res = 1 + ( fprime(w)/f(w) ) * ( OneMinusCapitalF(w) / f(w) ) 
  return(-res)
}


FromNumerical = function(w){
  OneMinusCapitalF(w) / f(w)  
}


# The derivative of ( (1-F(w))/f(w) ) w.r.t w at point w = 1 according to the theoretical derivation is:
w = seq(1, 5, by=1)
FromTheory(w)

# Numerical check for derivative of ( (1-F(w))/f(w) ) w.r.t w at point w = 1 is:
grad(func = FromNumerical, x = w)


abs( FromTheory(w) - grad(func = FromNumerical, x = w) )


