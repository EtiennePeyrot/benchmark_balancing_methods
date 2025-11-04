# library import

library(mvtnorm)

# function


GenerateData = function(n,
                        prop_trt = c("low", "moderate", "high"),
                        confdg_lvl = c("low", "moderate", "high")) {
  # Generate a samples of size n with various level of confounding and treatment
  # rarity.
  #
  # args :
  #       - n : integer, the sample size
  #       - prop_trt : str, the rarity of the treatement, can be
  #                    "low" (25%), "moderate" (50%), or "high" (75%).
  #       - confdg : str, the level of confounding for the simulation. can be
  #                  "low", "moderate" or "high".
  
  if (n <= 0 && n != floor(n))
    stop("'n' must be a positive integer")
  
  prop_trt = match.arg(prop_trt)
  confdg_lvl = match.arg(confdg_lvl)
  
  expit = function(f) 1 / (1 + exp(-f))
  
  # covariates
  
  Sigma = diag(rep(1,10)) # covariance matrix
  Sigma[5,1] = Sigma[1,5] = 0.2
  Sigma[8,3] = Sigma[3,8] = 0.2
  Sigma[6,2] = Sigma[2,6] = 0.9
  Sigma[9,4] = Sigma[4,9] = 0.9
  
  X = rmvnorm(n, rep(0, 10), Sigma)
  X[ ,c(1,3,5,6,8,9)] = as.numeric(X[ ,c(1,3,5,6,8,9)] > 0)
  
  # treatment assignment
  g = switch(confdg_lvl, "low" = 0.7, "moderate" = 1.2, "high" = 2.15)
  b0 = array(data = c(-1.48, -1.87, -2.73, -0.19, -0.32, -0.54, 1.07, 1.20, 1.59),
             dim = c(3,3),
             dimnames = list(confdf_lvl=c("low", "moderate", "high"),
                             prop_trt=c("low", "moderate", "high")))
  b0 = b0[confdg_lvl, prop_trt]
  b = c(0.8, -0.25, 0.6, -0.4, -0.8, -0.5, 0.7, 0, 0, 0)
  A = rbinom(n, 1, expit( b0 + g * (X %*% b + 0.5*X[,1]*X[,2]^2) ))
  
  
  # Outcome
  a0 = switch(confdg_lvl, "low" = -1.30, "moderate" = -1.65, "high" = -2.48)
  a = c(0.9, -1.08, -2.19, -0.6, 0, 0, 0, 0.71, -0.19, 0.26)
  Y = rbinom(n, 1, expit(a0 + g * (X %*% a + 0.5*X[,3]*X[,4]^2) ))
  
  
  # standardize continuous variable
  X = data.frame(X)
  X[,-c(1,3,5,6,8,9)] = scale(X[,-c(1,3,5,6,8,9)])
  
  output = list("X" = X,
                "A" = A,
                "Y" = Y,
                "confdg_lvl" = confdg_lvl,
                "prop_trt" = prop_trt)
  
  return(output)
}
