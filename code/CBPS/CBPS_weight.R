# import library

library(WeightIt)
if(!"rootSolve" %in% .packages(TRUE)) stop("package 'rootSolve' not installed")
library(rootSolve)

# functions


CBPS_weight = function(X, A, over = FALSE, twostep = NULL) {
  
  A = as.logical(A)
  
  fit = list(
    ATE = weightit(A ~ X, method = "cbps", estimand = "ATE", over = over,
                   twostep = twostep, include.obj = TRUE),
    ATT = weightit(A ~ X, method = "cbps", estimand = "ATT", over = over,
                   twostep = twostep, include.obj = TRUE, focal = TRUE)
  )
  
  W = lapply(fit, `[[`, "weights")
  
  # checks
  # W = lapply(fit, `[[`, "weights")
  # ps = lapply(fit, `[[`, "ps")
  # eta = lapply(fit, function(f) f$obj$par)
  # U = scale(svd(X)$u)
  # round(log10(abs(ps$ATE - plogis(cbind(1,U)%*% eta$ATE)))) |> table() |> print()
  # round(log10(abs(ps$ATT - plogis(cbind(1,U)%*% eta$ATT)))) |> table() |> print()
  # all.equal(W$ATE, ifelse(A, 1 / ps$ATE, 1 / (1 - ps$ATE))) |> print()
  # all.equal(W$ATT, ifelse(A, 1, ps$ATT / (1 - ps$ATT))) |> print()
  # browser()
  
  list(
    W = list(ATE = W$ATE / nrow(X), ATT = W$ATT / sum(A)),
    ps = lapply(fit, `[[`, "ps"),
    eta = lapply(fit, function(f) f$obj$par),
    U = scale(svd(X)$u),
    over = over
  )
}

CBPS_JI_weight = function(X, A) {
  CBPS_weight(X = X, A = A, over = FALSE)
}

CBPS_OI_weight = function(X, A, twostep = FALSE) {
  CBPS_weight(X = X, A = A, over = TRUE, twostep = twostep)
}
