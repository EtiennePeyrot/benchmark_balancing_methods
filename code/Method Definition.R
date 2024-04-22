# This script implement different balancing techniques and some nuisance
# functions which are going to be compared using the function "compute_method".

#### Initialization ####

debug_method_definition = FALSE
if (debug_method_definition) {
  setwd("G:/Mon Drive/Stage M2/Code")
  rm(list = ls()[ls() != "debug_method_definition"])
  graphics.off()
  set.seed(123)
}


# Library import

library(magrittr) # pipeline
library(randomForest) # random forest
library(caret) # regression + classification
library(WeightIt) # Energy balancing
library(osqp) # quadratic problem solver
library(sandwich) # var of coef in linear model
library(Mestim) # get var with M estim theory


##### Generic Functions #####

expit = function(x) 1/(1+exp(-x))

EstATE = function(Y, Z, W, y0, y1) {
  # Compute an estimation of the ATE for a given set of weight
  # args :
  #       - Y : vector containing the outcome
  #       - Z : logical vector containning treatment value
  #       - W: vector containing the weights
  #       - y0 / y1 : vector containing an estimation of the response surfaces
  
  ATE = mean(y1 - y0) + sum( W*ifelse(Z,(Y-y1), -(Y-y0)) )
  return(c("ATE" = sum((2*Z-1)*W*Y), "ATE.dr" = ATE))
}

EstATT = function(Y, Z, W, y0) {
  # Compute an estimation of the ATT for a given set of weight
  # args :
  #       - Y : vector containing the outcome
  #       - Z : logical vector containing treatment value
  #       - W : vector containing the weights
  #       - y0 : vector containing an estimation of the response surfaces
  if (!is.logical(Z)) Z = Z == 1
  ATT = mean(Y[Z]-y0[Z]) - sum(W*(Y[!Z]-y0[!Z]))
  return(c("ATT" = mean(Y[Z]) - sum(W*Y[!Z]), "ATT.dr" = ATT))
}

var_1nn = function(W, X, Y, Z) {
  if (!is.logical(Z))
    stop("'Z' lmust be logical.")
  sapply(list(Z, !Z), function(idx) {
    dist(X[idx, ], diag = T, upper = T) %>%
      as.matrix %>%
      `+`(diag(max(.) + 1, nrow = nrow(.), ncol = ncol(.))) %>%
      apply(1, which.min) %>%
      { .5 * (Y[idx] - Y[idx][.])^2 } %>%
      { sum(. * W[idx]^2) / sum(idx)^2 }
  }) %>% sum
}


standard_output = function(W.ATE, W.ATT, X, Y, y0, y1, Z, var,
                           mask.ATE = rep(T, length(Y)),
                           mask.ATT = rep(T, length(Y))) {
  if (!all(sapply(list(Z, mask.ATE, mask.ATT), is.logical)))
    stop(paste(paste(
      c("'Z'", "'mask.ATE'", "'mask.ATT'")[sapply(list(Z, mask.ATE, mask.ATT), is.logical)],
      collapse = ", "), "must be logical"))
  
  # compute the output returned by the method aka. AIPW, EB, KOM, TLF
  lm.ATE = lm(Y[mask.ATE] ~ Z[mask.ATE], weights = W.ATE[mask.ATE])
  lm.ATT = lm(Y[mask.ATT] ~ Z[mask.ATT], weights = W.ATT[mask.ATT])
  
  return( c(EstATE(Y = Y[mask.ATE], Z = Z[mask.ATE], W = W.ATE[mask.ATE], y0 = y0[mask.ATE], y1 = y1[mask.ATE]),
            "ATE.var" = if (var) var_1nn(W = W.ATE[mask.ATE], X = X[mask.ATE, ], Y = Y[mask.ATE], Z = Z[mask.ATE]) else NULL,
            "ATE.alt" = lm.ATE$coefficients[[2]],
            "ATE.alt.var" = if (var) diag(sandwich(lm.ATE))[[2]] else NULL,
            EstATT(Y = Y[mask.ATT], Z = Z[mask.ATT], W = W.ATT[mask.ATT & !Z], y0 = y0[mask.ATT]),
            "ATT.var" = if (var) var_1nn(W = W.ATT[mask.ATT], X = X[mask.ATT, ], Y = Y[mask.ATT], Z = Z[mask.ATT]) else NULL,
            "ATT.alt" = lm.ATT$coefficients[[2]],
            "ATT.alt.var" = if (var) diag(sandwich(lm.ATT))[[2]] else NULL) )
}

#### Nuisance Functions ####

Propensity_Score = function(X, Z, method, target = NULL,
                            trained.model = NULL, args = list()) {
  # Estimate the propensity score and return either the model or some prediction
  # args :
  #       - X : data.table or matrix  used to return a trained model or some prediction
  #       - Z : logical vector containing treatment value
  #       - method : char used to select the regression method to use
  #       - folds : integer setting the number of folds for the cross validation
  if (is.null(target)) target = "both"
  
  if (method == "rf") {
    # random forest model
    
    if (target != "predict") {
      # train the model on the data
      trt_ = factor(ifelse(Z, "test", "ctrl"), levels = c("ctrl", "test"))
      trctrl = trainControl(method="cv", number=5, classProbs=TRUE,
                            summaryFunction=mnLogLoss, allowParallel=FALSE)
      
      trained.model = train(
        trt_ ~ .,
        data = data.frame(X, "trt_" = trt_),
        method = "rf",
        metric = "logLoss",
        trControl = trctrl)
    }
    if (target == "train") return(trained.model)
    
    # predict the value of interest on the data
    return(list(pred = predict(trained.model, newdata = data.frame(X), type = "prob")[,"test"],
                model.name = "rf"))
  }
  
  
  
  if (method == "true") {
    # logistic regression with well specified model
    
    if (target != "predict") {
      # train the model on the data
      trained.model = glm(Z ~ ., data = cbind(X,X[,1]*X[,2]^2,Z), family = "binomial")
    }
    if (target == "train") return(trained.model)
    
    # predict the value of interest on the data
    return(list(
      "pred" = predict(trained.model, newdata = cbind(X,X[,1]*X[,2]^2), type = "response"),
      "model.coef" = `names<-`(coef(trained.model), NULL),
      "model.formula" = expression(cbind(1, X, X[ ,1] * X[ ,2]^2)),
      "model.name" = "log reg"))
  }
  
  
  
  if (method == "misspecified") {
    # logistic regression with well misspecified model
    
    if (target != "predict") {
      # train the model on the data
      trained.model = glm(Z ~ ., data = cbind(X,Z), family = "binomial")
    }
    if (target == "train") return(trained.model)
    
    # predict the value of interset on the data
    return(list(
      "pred" = predict(trained.model, newdata = X, type = "response"),
      "model.coef" = `names<-`(coef(trained.model), NULL),
      "model.formula" = expression(cbind(1,X)),
      "model.name" = "log reg"))
  }
}


Response_Surface = function(X, Y, Z, method, target = NULL,
                            trained.model = NULL, args = list()) {
  # Estimate the response surfaces
  # args :
  #       - data.train : data.table used as the training dataset
  #       - data.test : data.table used as the testing dataset
  #       - Z : logical vector containning treatment value
  #       - Y : vector containing the outcome
  #       - method : char used to select the regression method to use
  #       - folds : integer setting the number of folds for the cross validation
  if (is.null(target)) target = "both"
  
  if (method == "rf") {
    # random forest model
    
    if (target != "predict") {
      # train the model on the data
      
      y_ = factor(ifelse(Y, "alive", "dead"), levels = c("alive", "dead"))
      trctrl = trainControl(method="cv", number=5, classProbs=TRUE,
                            summaryFunction=mnLogLoss, allowParallel=FALSE)

      trained.model = list(
        "y0" = train(
          y_ ~ .,
          data = data.frame(X, "y_" = y_),
          subset = Z == 0,
          method = "rf",
          metric = "logLoss",
          trControl = trctrl),
        
        "y1" = train(
          y_ ~ .,
          data = data.frame(X, "y_" = y_),
          subset = Z == 1,
          method = "rf",
          metric = "logLoss",
          trControl = trctrl))
    }
    if (target == "train") return(trained.model)
    
    # predict the value of interest on the data
    return(lapply(trained.model, function(m) predict(m, newdata = X, type = "prob")[,"alive"]))
  }
  
  
  if (method == "true") {
    # logistic regression with well specified model
    
    if (target != "predict") {
    # train the model on the data
      trained.model = list(
        "y0" = glm(Y ~ ., data = cbind(X, X[, 3] * X[, 4] ^ 2, Y), subset = Z == 0, family = "binomial"),
        "y1" = glm(Y ~ ., data = cbind(X, X[, 3] * X[, 4] ^ 2, Y), subset = Z == 1, family = "binomial"))
    }
    if (target == "train") return(trained.model)
    
    # predict the value of interest on the data
    return(c(
      lapply(trained.model, predict, newdata = cbind(X, X[ ,3] * X[ ,4]^2), type = "response"),
      "model.coef" = lapply(trained.model, function(m) `names<-`(coef(m), NULL)),
      "model.formula" = list(expression(cbind(1, X, X[, 3] * X[, 4] ^ 2, Y))),
      "model.name" = "log reg"))
  }
  
  
  if (method == "misspecified") {
    # logistic regression with well misspecified model
    
    if (target != "predict") {
      # train the model on the data
      trained.model = list("y0" = glm(Y ~ ., data = cbind(X, Y), subset = Z == 0, family = "binomial"),
                           "y1" = glm(Y ~ ., data = cbind(X, Y), subset = Z == 1, family = "binomial"))
    }
    if (target == "train") return(trained.model)
    
    # predict the value of interset on the data
    return(c(
      lapply(trained.model, predict, newdata = X, type = "response"),
      "model.coef" = lapply(trained.model, function(m) `names<-`(coef(m), NULL)),
      "model.formula" = list(expression(cbind(1, X, Y))),
      "model.name" = "log reg"))
  }
}



#### AIPW ####

AIPW = function(X, Y, Z, propensity_score, response_surface, var = FALSE) {
  # Compute an estimation of ATE and ATT with AIPW
  # args :
  #       - X : X.frame containing the covariate
  #       - Y : vector containing the outcome
  #       - Z : logical vector containing treatment value
  #       - propensity_score : vector containing an estimation of the propensity score
  #       - response_surface : list of vector containing an estimation of the response surfaces
  ps =  propensity_score$pred
  if (!is.logical(Z)) Z = Z == 1
  y0 = response_surface$y0
  y1 = response_surface$y1
  
  # compute the weights
  W.ATE = 1 / ifelse(Z, ps, 1 - ps) / length(Z)
  W.ATT = ifelse(Z, 1, ps / (1 - ps)) / sum(Z)
  
  # remove extreme weights
  idx.ATE = W.ATE <= quantile(W.ATE, probs = .99, type = 7)
  idx.ATT = W.ATT <= quantile(W.ATT, probs = .99, type = 7)
  
  return(c(
    standard_output(W.ATE = W.ATE, W.ATT = W.ATT, X = X, Y = Y, y0 = y0, y1 = y1,
                    Z = Z, var = var, mask.ATE = idx.ATE, mask.ATT = idx.ATT)))
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
  
  return(standard_output(W.ATE = W.ATE, W.ATT = W.ATT, X = X, Y = Y,
                         y0 = y0, y1 = y1, Z = Z, var = var))
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

  return(standard_output(W.ATE = W.ATE, W.ATT = W.ATT, X = X, Y = Y,
                         y0 = y0, y1 = y1, Z = Z, var = var))
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
  
  if (is.null(gamma) | is.null(lambda)) {
    hyper_pram = TLF_hyper_param(X = X, Y = Y, Z = Z, nbr_fold = nbr_fold)
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

    
  return(standard_output(W.ATE = W.ATE, W.ATT = W.ATT, X = X, Y = Y,
                         y0 = y0, y1 = y1, Z = Z, var = var))
}

#### Test Functions ####


if (debug_method_definition) {
  sample_size = 1000
  X = data.frame(x1=1:sample_size,
                 x2=as.numeric(runif(sample_size)<0.5),
                 x3=runif(sample_size),
                 x4=rnorm(sample_size))
  Y = as.numeric(runif(sample_size) < 0.25)
  Z = as.numeric(runif(sample_size) < 0.75)
  split = floor(0.5 * sample_size)
  
  data.train = X[1:split,]
  data.test  = X[(split+1):sample_size, ]
  
  output = list()
  
  for (nui_model in c("true", "misspecified", "rf")) {
    ## Predictions
    
    e_pred = Propensity_Score(
      X = data.train,
      Z = Z[1:split],
      method = nui_model,
      target = "train") %>%
      
      Propensity_Score(
        X = data.test,
        method = nui_model,
        target = "predict",
        trained.model = .)
    
    pred = Response_Surface(
      X = data.train,
      Z = Z[1:split],
      Y = Y[1:split],
      method = nui_model,
      target = "train") %>%
      
      Response_Surface(
        X = data.test,
        Z = Z[-1:-split],
        Y = NULL,
        method = nui_model,
        target = "predict",
        trained.model = .)
    
    args = list("X" = data.test,
                "Z" = Z[-1:-split],
                "Y" = Y[-1:-split],
                "propensity_score" = e_pred,
                "response_surface" = pred,
                "gamma" = 100,
                "lambda" = 1e-3,
                "var" = TRUE)
    
    for (method in c("AIPW", "KOM", "EnergyB", "TLF"))
    {
      method_fct = switch(method,
                          "AIPW" = AIPW,
                          "KOM" = KOM,
                          "EnergyB" = EnergyB,
                          "TLF" = Tailored_loss_function)
      output[[paste(method, nui_model,sep = "_")]] =
        try(do.call(method_fct, args[names(args) %in% formalArgs(method_fct)]))
    }
  }
} else rm(debug_method_definition)
