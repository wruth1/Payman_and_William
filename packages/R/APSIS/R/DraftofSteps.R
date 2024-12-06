#
## f target density
#

#
# input a target density function,
# argument: what family we need to choose
# n: sample size
# iter: number of iteration
# stop rule criteria: some tolerance
# iter and stop: goes with some default values
#
# output ---> sample + proposal function
#

# List of options for family of proposals: Normal scale + Gamma


#
# Start with a target density function from user
# We have proposal density function N(\mu,1)
#


#
# We could have a function that can apply all families at the same time:
# Normal - Gamma - Cauchy
#


#
# Add diagnosis plots to look at the iteration vs value of parameters
#



# Zoom in one of the examples

#
# f: my target density function --- pass it to a function
# set the proposal family
# Run the optimization
#
# Function arguments:
# f: target density
# proposal_family
# sample size: n
#
# control = list() contains all the technical parts of the algorithm:
#   stopping rule,
#   step size,
#   tolerance,
#   max.iter
#   max.time
#

optimise_proposal = function(f, n, proposal_family = 'Normal', control = list() ){

  # Check if it is a function

  # Check the family proposal

  # Call the optimization

  # optimise_proposal_normal()
  # optimise_proposal_gamma()

  # Find the optimum value for parameters

  # Return:
  # Generated sample + Optimum Param Value + Proposal density ( g(x) ) + Final fitted Pareto tail

}
















