########################################################################
#
# Implements Kernel Optimal Matching (KOM) weights with a Gaussian (RBF)
# kernel for the ATE and ATT using the OSQP quadratic-programming solver.
#
# References:
#   - Kallus & al., “Optimal Estimation of Generalized Average Treatment
#     Effects using Kernel Optimal Matching (KOM).”
#   - Kallus, “Generalized Optimal Matching Methods for Causal Inference (GOM).”
#
########################################################################


# import library

library(osqp)

# function

Kgauss = function(X, scale){
  D = as.matrix(dist(X, upper = TRUE, diag = TRUE))^2
  exp(- scale^2 * D)
}

cross_Kgauss = function(X0, X1, scale) {
  # squared distances via ||u-v||^2 = ||u||^2 + ||v||^2 - 2 u·v
  X0sq = rowSums(X0^2); X1sq = rowSums(X1^2)
  D2 = outer(X0sq, X1sq, "+") - 2 * (X0 %*% t(X1))
  exp(- scale^2 * D2)
}

KOM_weight = function(X, A, lambda, scale) {
  A  = as.logical(A)
  X = as.matrix(X)
  n  = nrow(X); n1 = sum(A); n0 = n - n1
  
  # tight condition in solver to ensure weights are non negative
  osqp_setting = osqp::osqpSettings(
    verbose   = FALSE,
    polish    = TRUE,
    eps_abs   = 1e-8,
    eps_rel   = 1e-8,
    eps_prim_inf = 1e-12,
    eps_dual_inf = 1e-12,
    max_iter  = 1e6
  )
  
  out = list(ATE = numeric(n), ATT = rep(1 / n1, n))
  
  ## ATE
  # treated
  K = cross_Kgauss(X, X[A, ], scale = scale$trt )
  
  out$ATE[A] = osqp::solve_osqp(
    P = 2 * (K[A, ] + (lambda$trt + 1e-8) * diag(n1)), # ensure P is PSD
    q = (-2 / n) * as.numeric(colSums(K)),
    A = rbind(rep(1, n1), diag(n1)),
    l = c(1, rep(0, n1)),
    u = c(1, rep(1, n1)),
    pars = osqp_setting
  )$x
  
  # controls
  K = cross_Kgauss(X, X[!A, ], scale = scale$ctrl )
  
  out$ATE[!A] = osqp::solve_osqp(
    P = 2 * (K[!A, ] + (lambda$ctrl + 1e-8) * diag(n0)), # ensure P is PSD
    q = (-2 / n) * as.numeric(colSums(K)),
    A = rbind(rep(1, n0), diag(n0)),
    l = c(1, rep(0, n0)),
    u = c(1, rep(1, n0)),
    pars = osqp_setting
  )$x
  
  ## ATT
  # controls
  out$ATT[!A] = osqp::solve_osqp(
    P = 2 * (K[!A, ] + (lambda$ctrl + 1e-8) * diag(n0)),
    q = (-2 / n1) * as.numeric(colSums(K[A, ])),
    A = rbind(rep(1, n0), diag(n0)),
    l = c(1, rep(0, n0)),
    u = c(1, rep(1, n0)),
    pars = osqp_setting
  )$x
  
  return(out)
}