########################################################################
#
# Functions to compute normalized IPTW weights (ATE and ATT).
# Three models are used to estimate propensity scores: logistic regression with
# oracle term, plain logistic regression, and random forest.
#
########################################################################

# function


IPTW_logreg_weight = function(X, A) {
  n = nrow(X)
  
  # learn with logistic regression
  trained.model = glm(A ~ ., data = cbind(X, A), family = "binomial")
  
  # propensity score
  ps = as.vector(predict(trained.model, newdata = X, type = "response"))
  
  # weight
  weight = list(
    ATE = 1 / nrow(X) / ifelse(A, ps, 1 - ps),
    ATT = ifelse(A, 1, ps / (1 - ps)) / sum(A)
  )
  return(list(W = weight, model = trained.model))
}