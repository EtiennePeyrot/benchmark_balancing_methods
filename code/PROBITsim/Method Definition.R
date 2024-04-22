# This script implement different balancing techniques and some nuisance
# functions which are going to be compared using the function "compute_method".

# Library import

library(magrittr) # pipeline
library(randomForest) # random forest
library(caret) # regression + classification
library(WeightIt) # Energy balancing
library(osqp) # quadratic problem solver
library(sandwich) # var of coef in linear model
library(Mestim) # get var with M estim theory


#### AIPW ####

AIPW = function(X, Y, Z, propensity_score, response_surface, var = FALSE) {
  # Compute an estimation of ATE and ATT with AIPW
  # args :
  #       - X : X.frame containing the covariate
  #       - Y : vector containing the outcome
  #       - Z : logical vector containing treatment value
  #       - propensity_score : vector containing an estimation of the propensity score
  #       - response_surface : list of vector containing an estimation of the response surfaces
  
  ps =  propensity_score
  if (!is.logical(Z)) Z = Z == 1
  y0 = response_surface$y0
  y1 = response_surface$y1
  
  # compute the weights
  W.ATE = 1 / ifelse(Z, ps, 1 - ps) / length(Z)
  W.ATT = ifelse(Z, 1, ps / (1 - ps)) / sum(Z)
  
  return(list("W.ATE" = W.ATE, "W.ATT" = W.ATT))
}


#### Kernel Optimal Matching ####

minus_log_likelihood = function(param,K,Y) {
  # Evaluate the minus log likelihood of the outcome
  # code copied from Nathan Kallus github
  # cf : https://github.com/CausalML/KOM-SATE/blob/master/Case-Study/spine_functions_3_overploly_d.R
  # args :
  #       - param : vector of length 2
  #       - K : Kernel matrix
  #       - Y : vector containing the outcome
  sn2 = param[2] ^ 2
  gamma2 = param[1] ^ 2
  sample_size = length(Y)
  
  K = gamma2 * K + sn2 * diag(sample_size)
  K_inv = try(solve(K, tol=1e-21)
  )
  
  if(all(class(K_inv) != "try-error")){
    z = determinant(K, logarithm=TRUE)
    K_log_det = as.numeric((z$sign*z$modulus)) # log-determinant of K
    out = t(Y) %*% K_inv %*% Y + K_log_det
  }
  return(out)
}


KOMAutomaticSelectionLambda = function(Y,K) {
  # Compute an 'optimal' lambda according to Kallus papers
  # code copied from Nathan Kallus github
  # cf : https://github.com/CausalML/KOM-SATE/blob/master/Case-Study/spine_functions_3_overploly_d.R
  # args :
  #       - Y : vector containing the outcome
  #       - K : kernel matrix
  
  if (length(Y)<=1) return(0) # trivial case
  
  Y_centered = Y - mean(Y)
  res = try(optim(par=c(1, 1),
                  fn=minus_log_likelihood,
                  method=c("L-BFGS-B"),
                  lower=rep(1e-8,2),
                  hessian=TRUE,
                  control=list(trace=0, maxit=1000),
                  K=K,
                  Y=Y_centered) )
  
  out = ifelse(class(res)=="try-error", 1, (res$par[2]/res$par[1])^2) # if the optimization fails
  return(out)
}


GaussianKernel = function(X, sigma = 1) {
  # Compute the gaussian kernel matrix for a population X
  # args :
  #       - X : data.frame
  exp(-.5 * as.matrix(dist(x = X ,method = "euclidean", diag = T, upper = T))^2 / sigma)
}


KOM = function(X, Y, Z, response_surface, kernel = GaussianKernel,
               lbda=NA, var = FALSE) {
  # Compute an estimation of the ATE and ATT using Kernel Optimal Matching method created by Nathan Kallus
  # args :
  #       - X : X.frame containing the covariate
  #       - Z : vectors of boolean, 1 if the patient is treated 0 else
  #       - Y : vectors of double, the score/value of interst of the patient
  #       - response_surface : list of vector, an estimation of the surface response
  #       - kernel: function computing the kernel for this method
  #       - lbda : vector of doubles or NA for automatic selection, an hyperparameter
  
  if (!is.logical(Z)) Z = Z==1
  y0 = response_surface$y0
  y1 = response_surface$y1
  
  sample_size = length(Z)
  nbr_trt = sum(Z)
  nbr_ctrl = sample_size - nbr_trt
  
  # compute the kernel matrix
  K = kernel(X)
  
  if (length(lbda)==1) lbda = c(lbda,lbda)
  
  if (any(is.na(lbda))) { 
    # automatic selection of lambda as described in Kallus paper
    if (is.na(lbda[1])) lbda[1] = KOMAutomaticSelectionLambda(Y=Y[!Z],K=K[!Z, !Z])
    if (is.na(lbda[2])) lbda[2] = KOMAutomaticSelectionLambda(Y=Y[Z],K=K[Z, Z])
  }
  lbda0 = lbda[1]+1e-8
  lbda1 = lbda[2]+1e-8
  
  # compute the weights
  
  ## weights for T0 in ATT
  
  W.ATT = solve_osqp(
    P = 2 * (K[!Z, !Z] + lbda0 * diag(nbr_ctrl)),
    q = t(rep(-2 / nbr_trt, nbr_trt)) %*% K[Z, !Z],
    A = rbind(rep(1, nbr_ctrl), diag(nbr_ctrl)),
    l = c(1, rep(0, nbr_ctrl)),
    u = c(1, rep(Inf, nbr_ctrl)),
    pars = osqpSettings(verbose = FALSE))$x %>%
    sapply(max, 0)
  
  ## weights for T0 in ATE
  
  W0.ATE = solve_osqp(
    P = 2 * (K[!Z, !Z] + lbda0 * diag(nbr_ctrl)),
    q = t(rep(-2 / sample_size, sample_size)) %*% K[ , !Z],
    A = rbind(rep(1, nbr_ctrl), diag(nbr_ctrl)),
    l = c(1, rep(0, nbr_ctrl)),
    u = c(1, rep(Inf, nbr_ctrl)),
    pars = osqpSettings(verbose = FALSE))$x %>%
    sapply(max, 0)
  
  
  ## weights for T1 in ATE
  
  W1.ATE = solve_osqp(
    P = 2 * (K[Z, Z] + lbda1 * diag(nbr_trt)),
    q = t(rep(-2 / sample_size, sample_size)) %*% K[ , Z],
    A = rbind(rep(1, nbr_trt), diag(nbr_trt)),
    l = c(1, rep(0, nbr_trt)),
    u = c(1, rep(Inf, nbr_trt)),
    pars = osqpSettings(verbose = FALSE))$x %>%
    sapply(max, 0)
  
  W.ATE = rep(NA,sample_size); W.ATE[Z] = W1.ATE; W.ATE[!Z] = W0.ATE
  W.ATT = rep(1/nbr_trt, sample_size) %>% {.[!Z] = W.ATT; .}
  return(list("W.ATE" = W.ATE, "W.ATT" = W.ATT))
}


#### Energy Balancing ####

EnergyB = function(X, Y, Z, response_surface, var = FALSE) {
  # Compute the estimation of ATE and ATT with KOM
  # args :
  #       - X : data.frame containing the covariate
  #       - Y : vectors of double, the score/value of interst of the patient
  #       - Z : vectors of boolean, 1 if the patient is treated 0 else
  #       - response_surface : list of vector, an estimation of the surface response
  
  Z = Z==1
  y0 = response_surface$y0
  y1 = response_surface$y1
  
  # compute the weights
  W.ATE = weightit(Z ~ X,
                   method = "energy",
                   estimand = "ATE")$weights / ifelse(Z, sum(Z), sum(!Z)) # normalize the weights
  
  W.ATT = weightit(Z ~ X,
                   method = "energy",
                   estimand = "ATT",
                   focal=TRUE)$weights / sum(!Z) # normalize the weights
  return(list("W.ATE" = W.ATE, "W.ATT" = W.ATT))
}


#### Tailored Loss Function ####

TLF_hyper_param = function(X, Y, Z, nbr_fold = 5, Gamma = NULL, Lambda = NULL) {
  Gamma =  if(is.null(Gamma)) 10^(-3:3) else c(Gamma)
  Lambda = if(is.null(Lambda)) 10^(-3:3) else c(Lambda)
  fold_idx = sort(rep(1:nbr_fold, length.out = nrow(X)))
  Z_shift =  2 * Z - 1
  
  # select hyper parameter by cross validation
  hyper_param = lapply(Gamma, function(gamma) {
    K = exp(-gamma * as.matrix(dist(X, "manhattan", T, T)))
    
    lapply(Lambda, function(lambda) {
      c("gamma"  = gamma,
        "lambda" = lambda,
        "norm"   = mean(
          sapply(seq(nbr_fold), function(fold) {
            idx_train = fold_idx != fold
            idx_valid = fold_idx == fold
            
            # train on idx_train
            K.tr = K[idx_train, idx_train]
            Z_shift.tr = Z_shift[idx_train]
            # ATE
            optm = optim(
              par = matrix(0, 1, sum(idx_train)),
              method = "Nelder-Mead",
              fn = function(theta) {
                f_i = theta %*% K.tr
                v = expit(f_i)
                mean(Z_shift.tr * f_i - exp(-Z_shift.tr * f_i)) - lambda * sqrt(v %*% K.tr %*% t(v))
              })
            theta_star = optm$par
            f_star = optm$value
            
            # compute average norm of the gradient of the tailored loss function on the validation sample
            eps = 1e-3
            K.vl = K[idx_train, idx_valid]
            K.vl.norm = K[idx_valid, idx_valid]
            Z_shift.vl = Z_shift[idx_valid]
            norm = mean(
              apply(sweep(eps * diag(length(theta_star)), 1, theta_star, "+"), 2,
                    function(theta) {
                      f_i = theta %*% K.vl
                      v = expit(f_i)
                      mean(Z_shift.vl * f_i - exp(-Z_shift.vl * f_i)) - lambda * sqrt(v %*% K.vl.norm %*% t(v)) - f_star
                    })^2)
          })
        ))
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
  
  idx_param = which.min(hyper_param[ , "norm"])
  return(hyper_param[idx_param, 1:2])
}

Tailored_loss_function = function(X, Y, Z, response_surface, gamma = NULL,
                                  lambda = NULL, nbr_fold = NULL, var = FALSE) {
  # Compute the estimation of ATE and ATT with KOM
  # args :
  #       - X : data.frame containing the covariate
  #       - Y : vectors of double, the score/value of interst of the patient
  #       - Z : vectors of boolean or double/int, equal to 1 if the patient is treated 0 else
  #       - response_surface : list of vector, an estimation of the surface response
  #       - gamma : double, hyper parameter for the Laplacian kernel
  #       - lambda : double, hyper parameter for the penalization in the optimization problem
  #       - nbr_fold : int, number
  
  if (is.null(nbr_fold) & (is.null(gamma) | is.null(lambda)))
    stop("if 'nbr_fold' is NULL then 'gamma' and 'lambda' must be specified")
  
  if (!is.null(nbr_fold)) if (nbr_fold > nrow(X))
    stop("'nbr_fold' must be smaller or equal to 'nrow(X)'")
  
  if (!is.null(gamma)) if (gamma <= 0)
    stop("'gamma' must be positive")
  
  if (!is.logical(Z)) Z = Z==1
  
  y0 = response_surface$y0
  y1 = response_surface$y1
  Z_shift =  2 * Z - 1
  
  if (is.null(gamma) || is.null(lambda)) {
    hyper_param = TLF_hyper_param(X = X, Y = Y, Z = Z, nbr_fold = nbr_fold)
    gamma  = hyper_param["gamma" ]
    lambda = hyper_param["lambda"]
  }
  
  K = exp(-gamma * as.matrix(dist(X, "manhattan", T, T)))
  
  # ATE
  f_star = optim(par = matrix(0, 1, nrow(X)),
                 method = "Nelder-Mead",
                 fn = function(theta) {
                   f_i = theta %*% K
                   v = expit(f_i)
                   mean(Z_shift * f_i - exp(-Z_shift * f_i)) - lambda * sqrt(v %*% K %*% t(v))
                 })$par %*% K
  p_star = expit(f_star)
  W.ATE = ifelse(Z, 1/p_star, 1/(1-p_star))
  W.ATE = W.ATE / ifelse(Z, sum(W.ATE[Z]), sum(W.ATE[!Z]))
  
  # ATT
  f_star = optim(par = matrix(0, 1, nrow(X)),
                 method = "Nelder-Mead",
                 fn = function(theta) {
                   f_i = theta %*% K
                   v = expit(f_i)
                   mean(ifelse(Z, f_i, -exp(f_i))) - lambda * sqrt(v %*% K %*% t(v))
                 })$par %*% K
  p_star = expit(f_star)
  W.ATT = ifelse(Z, 1, p_star / ( 1 - p_star))
  W.ATT = W.ATT / ifelse(Z, sum(W.ATT[Z]), sum(W.ATT[!Z]))
  return(list("W.ATE" = W.ATE, "W.ATT" = W.ATT))
}