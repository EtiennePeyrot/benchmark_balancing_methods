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
  
  clip = function(p) pmin(pmax(p,0),1)
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
  b0_ = b0[confdg_lvl, prop_trt]
  
  b = c(0.8, -0.25, 0.6, -0.4, -0.8, -0.5, 0.7, 0, 0, 0)
  A = rbinom(n, 1, expit( b0_ + g * (X %*% b + 0.5*X[,1]*X[,2]^2) ))
  
  
  # Outcome
  a0 = c(-1.299, -1.645, -2.485, -1.3, -1.639, -2.487, -1.296, -1.649, -2.461,
         -1.382, -1.767, -2.653, -1.479, -1.877, -2.814, -1.558, -1.976, -2.975,
         -1.469, -1.878, -2.814, -1.65, -2.104, -3.155, -1.856, -2.375, -3.551)
  a0 = array(data = a0, dim = rep(3,3),
             dimnames = c(dimnames(b0),list(effect=c("none","small","high"))))
  
  a = c(0.9, -1.08, -2.19, -0.6, 0, 0, 0, 0.71, -0.19, 0.26)
  y_lin = g * (X %*% a + 0.5*X[,3]*X[,4]^2)
  Y = lapply(c("none", "small", "high"),
             function(effect) {
               ate = switch(effect, "none"=0.00, "small"=0.05, "high"=0.10)
               a0_ = a0[confdg_lvl, prop_trt, effect]
               prob_event = clip(ate * A + expit((a0_ + y_lin)))
               rbinom(n, 1, prob_event)
             })
  names(Y) = c("none", "small", "high")
  
  # standardize continuous variable
  X = data.frame(X)
  X[ ,-c(1,3,5,6,8,9)] = scale(X[,-c(1,3,5,6,8,9)])
  
  output = list("X" = X,
                "A" = A,
                "Y" = Y,
                "confdg_lvl" = confdg_lvl,
                "prop_trt" = prop_trt)
  
  return(output)
}
