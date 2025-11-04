kernel.balance <- function(T, X, K = NULL, d = NULL,
                           alpha, beta, lambda,
                           r = ncol(X),
                           intercept = TRUE, # add this later
                           normalize = TRUE,
                           maxiter = 100, tol = 1e-4,
                           eta.init = NULL,
                           verbose = FALSE) {
  
  if (alpha > 0 || alpha < -1 || beta > 0 || beta < -1) {
    stop("alpha and beta must be between -1 and 0.")
  }
  
  if (!is.null(K)) {
    K.eigen <- eigen(K)
    U <- K.eigen$vectors
    d <- K.eigen$values
    if (is.null(r)) {
      r <- sum(d > 1e-12)
    }
    U <- U[, 1:r]
    d <- 1/d[1:r]
  } else {
    U <- X
    if (is.null(d)) {
      d <- rep(1, r)
    }
  }
  
  if (intercept) {
    U <- cbind(rep(1, length(T)), U)
    d <- c(0, d)
  }
  
  n <- length(T)
  
  S.func <- function(f) {
    prob <- link(f)
    S <- rep(0, length(f))
    if (alpha == -1 && beta == -1) {
      S[T == 1] <- f[T == 1] - 1 / prob[T == 1]
      S[T == 0] <- - f[T == 0] - 1 / (1 - prob[T == 0])
    } else if (alpha == -1 && beta == 0) {
      S[T == 1] <- - 1 / prob[T == 1]
      S[T == 0] <- - f[T == 0]
    } else if (alpha == 0 && beta == -1) {
      S[T == 1] <- f[T == 1]
      S[T == 0] <- - 1 / (1 - prob[T == 0])
    } else if (alpha == 0 && beta == 0) {
      S[T == 1] <- log(prob[T == 1])
      S[T == 0] <- log(1 - prob[T == 0])
    } else {
      for (i in 1:length(prob)) {
        S[i] <- integrate(function(p) (T[i] - p) * p^(alpha - 1) * (1 - p)^(beta - 1), lower = 1/2, upper = prob[i])$value
      }
    }
    sum(S) / n
  }
  
  S.prime <- function(f) {
    prob <- link(f)
    Sp <- rep(0, length(f))
    Sp[T == 1] <- prob[T == 1]^alpha * (1 - prob[T == 1])^(beta + 1)
    Sp[T == 0] <- - prob[T == 0]^(alpha + 1) * (1 - prob[T == 0])^beta
    Sp
  }
  
  S.prime.prime <- function(f) {
    p <- link(f)
    Spp <- rep(0, length(p))
    Spp[T == 1] <- (alpha * p[T == 1]^alpha * (1 - p[T == 1])^(beta + 2) - (beta + 1) * p[T == 1]^(alpha + 1) * (1 - p[T == 1])^(beta + 1))
    Spp[T == 0] <- (beta * (1 - p[T == 0])^beta * p[T == 0]^(alpha + 2) - (alpha + 1) * (1 - p[T == 0])^(beta + 1) * p[T == 0]^(alpha + 1))
    Spp
  }
  
  link <- function(f) {
    1 / (1 + exp(-f))
  }
  
  gradient <- function(eta) {
    f <- U %*% eta
    t(U) %*% S.prime(f) / n  - 2 * lambda * d * eta
  }
  
  hessian <- function(eta) {
    f <- U %*% eta
    t(U * as.vector(S.prime.prime(f))) %*% U / n - 2 * lambda * diag(d)
  }
  
  n <- length(T)
  ## initialize
  if (is.null(eta.init)) {
    eta.init <- rep(0, r + 1)
  }
  eta <- eta.init
  
  get.obj.value <- function(eta) {
    S.func(U %*% eta) - lambda * sum(eta[-1]^2 * d[-1])
  }
  
  ## K.inv <- solve(K)
  
  for (iter in 1:maxiter) {
    
    ## print(min(plogis(U %*% eta)))
    ## print(max(plogis(U %*% eta)))
    obj.value <- get.obj.value(eta)
    if (verbose) {
      print(paste("Objective:", obj.value))
    }
    
    grad <- gradient(eta)
    hess <- hessian(eta)
    ## print(kappa(hess))
    ## print(kappa(hess %*% K.inv))
    ## print(kappa(K.inv %*% hess %*% K.inv))
    library(MASS)
    ## hess.inv <- ginv(hess)
    ## newton <- hess.inv %*% grad
    newton <- solve(hess, grad)
    ## print(kappa(hess %*% solve(K)))
    
    if (verbose) {
      print(paste("Max gradient:", max(abs(grad))))
    }
    if (sum(grad^2) < tol) {
      break
    }
    
    eta.new <- eta - newton
    if (verbose) {
      print(paste("Newton value:", get.obj.value(eta.new)))
    }
    if (get.obj.value(eta.new) > obj.value) {
      if (verbose) {
        print("Newton step")
      }
      eta <- eta.new
    } else {
      m <- -Inf
      upper <- 1
      while (m == -Inf) {
        line.optimizer <- optimize(f = function(step.size) get.obj.value(eta - step.size * newton), c(0, upper), maximum = TRUE)
        m <- line.optimizer$objective
        upper <- upper / 5
      }
      if (line.optimizer$objective > obj.value) {
        if (verbose) {
          print("Partial Newton step")
        }
        eta <- eta - line.optimizer$maximum * newton
      } else {
        if (verbose) {
          print("Gradient step")
        }
        line.optimizer <- optimize(f = function(step.size) get.obj.value(eta + step.size * grad), c(0, 1), maximum = TRUE)
        if (verbose) {
          print(line.optimizer)
        }
        if (line.optimizer$objective > obj.value) {
          eta <- eta + line.optimizer$maximum * grad
        } else {
          break
        }
      }
    }
  }
  
  if (sum(grad^2) < tol) {
    converged <- TRUE
  } else {
    converged <- FALSE
  }
  
  f <- U %*% eta
  prob <- link(f)
  weights <- compute.weights(T, prob, alpha, beta, normalize)
  
  return(list(eta = eta,
              f = f,
              S = S.func(f),
              p = prob,
              grad = grad,
              converged = converged,
              weights = weights))
  
}

kernel.balance.path <- function(T, X, K = NULL, d = NULL,
                                alpha, beta, lambda,
                                r = ncol(X),
                                intercept = TRUE, # add this later
                                normalize = TRUE,
                                maxiter = 100, tol = 1e-4,
                                eta.init = NULL,
                                verbose = FALSE) {
  
  if (alpha > 0 || alpha < -1 || beta > 0 || beta < -1) {
    stop("alpha and beta must be between -1 and 0.")
  }
  
  if (!is.null(K)) {
    K.eigen <- eigen(K)
    U <- K.eigen$vectors
    d <- K.eigen$values
    if (is.null(r)) {
      r <- sum(d > 1e-12)
    }
    U <- U[, 1:r]
    d <- 1/d[1:r]
  } else {
    U <- X
    if (is.null(d)) {
      d <- rep(1, r)
    }
  }
  
  if (is.null(eta.init)) {
    eta.init <- rep(0, r + 1)
  }
  
  result <- list()
  for (i in 1:length(lambda)) {
    result[[i]] <- kernel.balance(T, U, d = d,
                                  alpha = alpha, beta = beta,
                                  lambda = lambda[i],
                                  intercept = intercept,
                                  normalize = normalize,
                                  maxiter = maxiter,
                                  tol = tol,
                                  eta.init = eta.init,
                                  verbose = verbose)
    eta.init <- result[[i]]$eta
  }
  
  result
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
