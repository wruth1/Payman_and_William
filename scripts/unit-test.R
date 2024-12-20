#
# This script tests values for f(w) and its derivative by comparing results from theory to 
# results from numerical computation.
# 

# Load libraries for Lambert W function and numerical derivations
library(pracma)
library(numDeriv)

# Change default digits accuracy
options(digits = 18)


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


lambertWp2prime = function(w){
  LLWp <- lambertWp(x=w)
  num <- 2 + LLWp
  den <- 1 + LLWp
  prime <-  LLWp / (w * (1 + LLWp) )
  res <- (-1) * prime * prime * (num/den)
  return(res)
}


lambertWn2prime = function(w){
  LLWn <- lambertWn(x=w)
  num <- 2 + LLWn
  den <- 1 + LLWn
  prime <-  LLWn / (w * (1 + LLWn) )
  res <- (-1) * prime * prime * (num/den)
  return(res)
}


lambertWpDeriv = function(w){
  LLWp <- lambertWp(x=w)
  num  <- LLWp
  den  <- w * (1 + LLWp)
  res  <- num / den
  return(res)
}


lambertWnDeriv = function(w){
  print(w)
  LLWn <- lambertWn(x=w)
  num  <- LLWn
  den  <- w * (1 + LLWn)
  res  <- num / den
  return(res)
}


# Defined in formula 40 in section 2.3 of pdf file
SurvivalFunction = function(w, lambda, gamma, alpha, c){
  
  res <- exp( -(Xlw(w, gamma, alpha, c)/lambda) ) - exp( -(Xuw(w, gamma, alpha, c)/lambda) )
  return(res)
  
}


# Defined in formula 54 in section 2.3 of pdf file
f = function(w, lambda, A, gamma, alpha, c){
  
  eval_point <- A * w^{1 / (alpha-1)}
  term1 <- ( lambertWn(eval_point) / (1 + lambertWn(eval_point)) ) * exp( -(Xuw(w, gamma, alpha, c)/lambda) )
  term2 <- ( lambertWp(eval_point) / (1 + lambertWp(eval_point)) ) * exp( -(Xlw(w, gamma, alpha, c)/lambda) )
  output = gamma * (term1 - term2) / (lambda * w)
  return(output)
  
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

# Defined in formula 79 of Gamma LRT.pdf
fprime = function(w, lambda, A, gamma, alpha, c){
  
  eval_point <- A * w^( 1 / (alpha-1) )
  
  mult_term <- (A/(alpha-1)) * w^( (2-alpha) / (alpha-1) )
  
  term1 <- (gamma/lambda) * exp( -(Xuw(w, gamma, alpha, c)/lambda) ) * lambertWn2prime(eval_point) * mult_term
  term2 <- (gamma^2/lambda^2) * exp( -(Xuw(w, gamma, alpha, c)/lambda) ) * lambertWnDeriv(eval_point)
  fprime_u <- term1 - term2

  term1 <- (gamma/lambda) * exp( -(Xlw(w, gamma, alpha, c)/lambda) ) * lambertWp2prime(eval_point) * mult_term
  term2 <- (gamma^2/lambda^2) * exp( -(Xlw(w, gamma, alpha, c)/lambda) ) * lambertWpDeriv(eval_point)
  fprime_l <- term1 - term2
  
  return(fprime_u - fprime_l)
}

# Defined in Gamma LRT Appendix.pdf
FromTheory = function(w, A, gamma, alpha, c){
  res = 1 + ( fprime(w, A, gamma, alpha, c)/f(w, A, gamma, alpha, c) ) * ( SurvivalFunction(w, gamma, alpha, c) / f(w, A, gamma, alpha, c) ) 
  return(-res)
}


FromNumerical = function(w, A, gamma, alpha, c){
  SurvivalFunction(w, gamma, alpha, c) / f(w, A, gamma, alpha, c)  
}



# Set some of the parameters to constant values
lambda <- 2
beta   <- 1
alpha  <- 3
gamma  <- 1/(1/beta - 1/lambda)
c      <- lambda / (beta^alpha * gamma(alpha))
A      <- (-1)/(gamma * (alpha-1) * c^(alpha-1))
w      <- 1
w_star <- c * (gamma * (alpha - 1) * exp(-1))^(alpha - 1)
some_Ws <- seq(0.1, 2, by = 0.01)


#
# Test 1: If f(w) equals to negative derivative of survival function
#

# We can check it with more parameter settings too

math_calc <- f(some_Ws, lambda, A, gamma, alpha, c)
num_grad  <- (-1) * grad(func = SurvivalFunction, x = some_Ws, lambda = lambda, gamma = gamma, alpha = alpha, c = c)

testthat::expect_equal(num_grad, math_calc)
math_calc - num_grad


#
# Test 2: If fprime(w) equals to derivative of f(w)
#

math_calc <- fprime(w = some_Ws, lambda = lambda, A = A, gamma = gamma, alpha = alpha, c = c)
num_grad  <- grad(func = f, x = some_Ws, lambda = lambda, A = A, gamma = gamma, alpha = alpha, c = c)

testthat::expect_equal(num_grad, math_calc)




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



