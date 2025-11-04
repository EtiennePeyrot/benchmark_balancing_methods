########################################################################
#
# Implements hyper-parameter selection for Kernel Optimal Matching (KOM)
# with a Gaussian (RBF) kernel for the ATE and ATT.
#
# Fit GP (Type-II ML) within each arm on (X_t, Y_t) -> (gamma^2, scale, sigma^2)
# Set lambda_t = sigma^2 / gamma^2; use only scale when building K.
#
# References:
#   - Kallus & al., “Optimal Estimation of Generalized Average Treatment
#     Effects using Kernel Optimal Matching (KOM).”
#   - Kallus, “Generalized Optimal Matching Methods for Causal Inference (GOM).”
#
########################################################################

# negative log likelihood
kom_nll = function(par_log, D2, y, grad=T) {
  n = length(y)
  
  gamma2 = exp(par_log[1])
  scale  = exp(par_log[2])
  sigma2 = exp(par_log[3])
  
  # Base Gaussian Gram without amplitude
  G  = exp(- (scale^2) * D2)
  Kp = gamma2 * G + sigma2 * diag(n)
  
  # Stable Cholesky
  R = try(chol(Kp), silent = TRUE)
  if (inherits(R, "try-error")) {
    eps = 1e-10 * mean(diag(Kp))
    R   = chol(Kp + eps * diag(n))
  }
  
  # α = K^{-1} y, log|K|
  alpha  = backsolve(R, forwardsolve(t(R), y))
  yKy    = sum(y * alpha)
  logdet = 2 * sum(log(diag(R)))
  nll    = yKy + logdet
  
  if(!grad) return(nll)
  
  # For gradients: use Ki = K^{-1} via Cholesky (full symmetric inverse)
  Ki = chol2inv(R)   # O(n^3)
  
  # Common sandwich: S = K^{-1} - αα^T
  # (∂/∂θ) [ y^T K^{-1} y + log|K| ] = tr( S * ∂K/∂θ )
  S  = Ki - tcrossprod(alpha)
  
  # Chain rule to log-parameters: d/d log θ = θ * d/d θ
  grad_log_gamma2 = gamma2 * sum(S * G)
  grad_log_scale  = -2 * gamma2 * scale^2  * sum(S * G * D2)
  grad_log_sigma2 = sigma2 * sum(diag(S))
  
  list(value = nll,
       grad  = c(grad_log_gamma2, grad_log_scale, grad_log_sigma2))
}

fit_GP = function(X, Y) {
  D2 = as.matrix(dist(X, diag = T, upper = T))^2
  
  fn = function(p)  kom_nll(p, D2, Y, grad=FALSE)
  gr = function(p)  kom_nll(p, D2, Y)$grad
  
  par = try(optim(par = c(gamma2 = 0, scale = 0, sigma2 = 0),
                  fn = fn, gr = gr,
                  method = "L-BFGS-B", lower = log(rep(1e-8, 3)),
                  control = list(trace = 0, maxit = 1000))$par,
            silent = FALSE)
  
  if (inherits(par, "try-error")) return(c(lambda = 1, scale = 1))
  
  par = exp(par) # remove log scale
  c(lambda = par[["sigma2"]] / par[["gamma2"]], scale = par[["scale"]])
}

KOM_hp = function(X, A, Y) {
  X = as.matrix(X)
  A = as.logical(A)
  
  list(
    ctrl = fit_GP(X[!A, ], Y[!A]),
    trt  = fit_GP(X[ A, ], Y[ A])
  )
}