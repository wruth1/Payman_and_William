
setwd("packages/R/APSIS/")

devtools::load_all()



target = function(x) dnorm(x)

n = 1000

mu0 = 1


optimise_proposal_normal(target, n, mu0)



update_proposal_normal(n, mu0, target, 1)

estimate_gradient_normal()



sum(weight)^2 / sum(weight^2)


x = xsample
w = weight
