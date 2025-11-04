kernel.balance = function(T, X, K = NULL, d = NULL,
                           alpha, beta, lambda,
                           r = NULL, K.eigen = NULL,
                           intercept = TRUE, # add this later
                           normalize = TRUE,
                           maxiter = 100, tol = 1e-4,
                           eta.init = NULL, energy = .9,
                           verbose = FALSE) {
  
  if (alpha > 0 || alpha < -1 || beta > 0 || beta < -1) {
    stop("alpha and beta must be between -1 and 0.")
  }
  
  n = length(T)
  if (is.null(K.eigen) & !is.null(K)) K.eigen = eigen(K,symmetric=TRUE)
  if (!is.null(K.eigen)) {
    d = K.eigen$values
    if (is.null(r)) r = max(ncol(X),which(cumsum(d/sum(d)) < energy))
    U = K.eigen$vectors[, 1:r]
    d = 1/d[1:r]
  }
  if (is.null(K.eigen) &  is.null(K)) {
    U = X
    if (is.null(d)) d = rep(1, r)
  }
  
  if (intercept) {
    U = cbind(1, U)
    d = c(0, d)
  }
  
  
  link = function(f) {
    1 / (1 + exp(-f))
  }
  
  S.func = function(f) {
    prob = link(f)
    S = rep(0, length(f))
    if (alpha == -1 && beta == -1) {
      S[T == 1] = f[T == 1] - 1 / prob[T == 1]
      S[T == 0] = - f[T == 0] - 1 / (1 - prob[T == 0])
    } else if (alpha == -1 && beta == 0) {
      S[T == 1] = - 1 / prob[T == 1]
      S[T == 0] = - f[T == 0]
    } else if (alpha == 0 && beta == -1) {
      S[T == 1] = f[T == 1]
      S[T == 0] = - 1 / (1 - prob[T == 0])
    } else if (alpha == 0 && beta == 0) {
      S[T == 1] = log(prob[T == 1])
      S[T == 0] = log(1 - prob[T == 0])
    } else {
      for (i in 1:length(prob)) {
        S[i] = integrate(function(p) (T[i] - p) * p^(alpha - 1) * (1 - p)^(beta - 1), lower = 1/2, upper = prob[i])$value
      }
    }
    sum(S) / n
  }
  
  S.prime = function(f) {
    prob = link(f)
    Sp = rep(0, length(f))
    Sp[T == 1] = prob[T == 1]^alpha * (1 - prob[T == 1])^(beta + 1)
    Sp[T == 0] = - prob[T == 0]^(alpha + 1) * (1 - prob[T == 0])^beta
    Sp
  }
  
  S.prime.prime = function(f) {
    p = link(f)
    Spp = rep(0, length(p))
    Spp[T == 1] = (alpha * p[T == 1]^alpha * (1 - p[T == 1])^(beta + 2) - (beta + 1) * p[T == 1]^(alpha + 1) * (1 - p[T == 1])^(beta + 1))
    Spp[T == 0] = (beta * (1 - p[T == 0])^beta * p[T == 0]^(alpha + 2) - (alpha + 1) * (1 - p[T == 0])^(beta + 1) * p[T == 0]^(alpha + 1))
    Spp
  }
  
  gradient = function(eta) {
    f = U %*% eta
    t(U) %*% S.prime(f) / n  - 2 * lambda * d * eta
  }
  
  hessian = function(eta) {
    f = U %*% eta
    t(U * as.vector(S.prime.prime(f))) %*% U / n - 2 * lambda * diag(d)
  }
  
  ## initialize
  if (is.null(eta.init)) eta.init = rep(0, r + intercept)
  eta = eta.init
  
  get.obj.value = function(eta) {
    S.func(U %*% eta) - lambda * sum(if(intercept) eta[-1]^2 * d[-1] else eta^2 * d)
  }
  
  
  if(verbose) cat("\n\noptimization start\n")
  converged = FALSE
  for (iter in 1:maxiter) {
    
    obj.value = get.obj.value(eta)
    if (verbose) cat("\n\nobjective:", obj.value)
    
    grad = gradient(eta)
    if (verbose) cat("\nmax gradient:", max(abs(grad)))
    if (sum(grad^2) < tol) {
      converged = TRUE
      break
    }
    hess = hessian(eta)
    do.gradiant = inherits(try(newton <- solve(hess, grad), silent=TRUE), "try-error")
    
    if(!do.gradiant) {
      # newton step
      obj.value.newton = get.obj.value(eta - newton)
      if (obj.value.newton > obj.value) {
        if (verbose) cat("\nnewton step")
        eta = eta - newton
        next
      }
      
      
      # partial newton step
      m = -Inf
      upper = 1
      while (m == -Inf) {
        line.optimizer = optimize(f = function(step.size) get.obj.value(eta - step.size * newton),
                                   interval = c(0, upper), maximum = TRUE)
        m = line.optimizer$objective
        upper = upper / 5
      }
      if (line.optimizer$objective > obj.value) {
        if (verbose) cat("\npartial Newton step")
        eta = eta - line.optimizer$maximum * newton
        next
      }
    }
    
    # gradient step
    line.optimizer = optimize(f = function(step.size) get.obj.value(eta + step.size * grad),
                               interval = c(0, 1), maximum = TRUE)
    
    
    
    if (line.optimizer$objective > obj.value) {
      if (verbose) cat("\ngradient step")
      eta = eta + line.optimizer$maximum * grad
      next
    }
    
    # fail
    break
  }
  if(verbose) cat("\n\noptimization end\n")
  
  f = U %*% eta
  prob = link(f)
  weights = compute.weights(T, prob, alpha, beta, normalize)
  
  proj = if (intercept) U[ ,-1] %*% diag(d[-1]) %*% eta[-1] else U %*% diag(d) %*% eta
  cst = if (intercept) eta[1] else 0
  predict = function(k_new) {
    link(unname(drop(k_new %*% proj  + cst)))
  }
  return(list(eta = eta,
              f = f,
              U = U,
              d = d,
              predict = predict,
              S = S.func(f),
              p = prob,
              grad = grad,
              converged = converged,
              weights = weights))
  
}

TLF_weight = function(X, A, sigma = NULL, lambda,
                      K.eigen = NULL, intercept = TRUE,
                      maxiter = 100, tol = 1e-4, energy = .9,
                      verbose = FALSE) {
  if (is.null(X) & is.null(sigma) & is.null(K.eigen)) stop("'X' and 'sigma' must be specified if 'K.eigen' is NULL")
  
  if (is.null(K.eigen)) {
    K = exp(-sigma * as.matrix(dist(x = X, method = "manhattan", diag = T, upper = T)))
    K.eigen = eigen(K)
  }
  
  res_ate = kernel.balance(
    T = A, X = X, K = NULL, K.eigen = K.eigen,
    alpha = -1, beta = -1,
    lambda = lambda$ATE,
    intercept = intercept, normalize = F,
    energy = .9, verbose = verbose
  )
  
  res_att = kernel.balance(
    T = A, X = X, K = NULL, K.eigen = K.eigen,
    alpha = 0, beta = -1,
    lambda = lambda$ATT,
    intercept = T, normalize = F,
    energy = .9, verbose = verbose
  )
  
  out = list(
    weights = list(ATE = res_ate$weights / length(A),
                   ATT = res_att$weights / sum(A)),
    model   = list(ATE = res_ate,
                   ATT = res_att)
  )
  
  return(out)
}

kernel.hp.gridsearch = function(T, X, alphas, betas,
                                lambdas, sigmas = NULL,
                                nbr_fold, energy = .9, intercept = TRUE,
                                maxiter = 100, tol = 1e-4,
                                verbose = FALSE, best = F, seed=NULL) {
  n = nrow(X)
  expit = function(x) 1 / (1 + exp(-x))
  
  if(length(alphas) != length(betas))
    stop("'alphas' and 'betas' length are not equal")
  
  D = as.matrix(dist(X, "manhattan", diag = T, upper = T))
  
  if(is.null(sigmas)) sigmas = 1/median(D[D>0])
  
  if(!is.null(seed)) set.seed(seed)
  fold = sample(rep(1:nbr_fold, length.out = n))
  
  S.prime = function(prob, T, alpha, beta) {
    Sp = rep(0, length(T))
    Sp[T == 1] = prob[T == 1]^alpha * (1 - prob[T == 1])^(beta + 1)
    Sp[T == 0] = - prob[T == 0]^(alpha + 1) * (1 - prob[T == 0])^beta
    Sp
  }
  
  
  pars = expand.grid(estimand_idx = 1:length(alphas), lambda = lambdas, sigma = sigmas)
  res = array(
    data = NA_real_,
    dim = c(length(lambdas), length(sigmas), length(alphas)),
    dimnames = list(lambda=lambdas, sigma=sigmas, estimand=paste0("alpha=",alphas,"; beta=",betas)))
  
  if(verbose) cat("\nTLF hyper param grid search for alpha =", alphas, ", beta =", betas, "\n")
  width = floor(log10(nrow(pars)))
  sigma_ = max(sigmas) + 1
  for (i in 1:nrow(pars)) {
    estimand_idx = pars$estimand_idx[i]; lambda = pars$lambda[i]; sigma = pars$sigma[i]
    alpha_beta_key = paste0("alpha=",alphas[estimand_idx],"; beta=",betas[estimand_idx])
    metric = rep(NA_real_, nbr_fold)
    
    if(verbose) cat(paste0("\riteration ",format(i,width=width),"/",nrow(pars)))
    
    # Compute Laplace kernel matrix K and eigen vectors if sigma changes
    if (sigma != sigma_) {
      K = exp(-sigma * D)
      sigma_ = sigma
      K.eigen = lapply(
        1:nbr_fold,
        function(idx_fold) eigen(K[fold != idx_fold, fold != idx_fold])
      )
    }
    
    for (idx_fold in 1:nbr_fold) {
      idx_tr  = fold != idx_fold; n_tr  = sum(idx_tr )
      idx_val = fold == idx_fold; n_val = sum(idx_val)
      
      # use train fold to predict the propensity score on validation fold
      opt = kernel.balance(
        T = T[idx_tr], X = X[idx_tr, ], K.eigen = K.eigen[[idx_fold]],
        alpha = alphas[estimand_idx], beta = betas[estimand_idx],
        lambda = lambda, tol = tol,
        intercept = intercept, normalize = F,
        energy = energy
      )
      
      
      if (!opt$converged) next
      # predict propensity score on validation fold
      p_val = opt$predict(K[idx_val, idx_tr])
      
      # compute unpenalized gradient
      dS = S.prime(prob = p_val, T = T[idx_val], alpha = alphas[estimand_idx], beta = betas[estimand_idx])
      
      metric[idx_fold] = mean(dS^2)
    }
    
    res[as.character(lambda), as.character(sigma), alpha_beta_key] = mean(metric)
  }
  if(verbose) cat("\nDone!\n")
  
  if (!best) return(drop(res))
  apply(res, 3, function(r) {
    par = which(r == min(r, na.rm=T), arr.ind = T)
    return(c("lambda" = lambdas[par[1]], "sigma" = sigmas[par[2]]))
  }, simplify = FALSE)
}

compute.weights <- function(T, prob, alpha, beta, normalize = TRUE) {
  
  weights <- rep(0, length(T))
  weights[T == 1] <- prob[T == 1]^alpha * (1 - prob[T == 1])^(beta + 1)
  weights[T == 0] <- prob[T == 0]^(alpha + 1) * (1 - prob[T == 0])^beta
  if (normalize) {
    n <- length(T)
    weights[T == 1] <- n * weights[T == 1] / sum(weights[T == 1])
    weights[T == 0] <- n * weights[T == 0] / sum(weights[T == 0])
  }
  weights
}

