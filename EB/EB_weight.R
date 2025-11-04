########################################################################
#
# Function to compute normalized Energy balancing weight weights (ATE and ATT).
# Rely on the implementation provided in the package WeightIt by Noah Greifer.
#
# Ref: Huling, J. D., & Mak, S. (2020). Energy Balancing of Covariate
# Distributions.
#
########################################################################


# import library

library(WeightIt)

# function

EB = function(X, A) {
  A = as.logical(A)
  
  norm = ifelse(A, 1 / sum(A), 1 / sum(1-A))
  # compute the weights
  weight = list(
    ATE = norm * weightit(A ~ X, method = "energy", estimand = "ATE")$weights,
    ATT = norm * weightit(A ~ X, method = "energy", estimand = "ATT", focal = T)$weights
  )
  
  return(weight)
}