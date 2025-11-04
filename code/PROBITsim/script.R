rm(list = ls())
graphics.off()
set.seed(123)

# file import

source("../EB/EB_weight.R")
source("../IPTW/IPTW_weight.R")
source("../IPTW/confidence_interval.R")
source("../KOM/KOM_weight.R")
source("../KOM/KOM_hp_selec.R")
source("../KOM/confidence_interval.R")
source("../TLF/TLF_weight.R")
source("../TLF/confidence_interval.R")
source("confidence_interval.R")

# library import

library(foreign) # import dta file

# data import

probitsim = as.data.frame(read.dta("PROBITsim.dta", convert.factors = FALSE))

#### Functions ####

.mycombine = function(l) {
  name = names(l[[1]])
  out = lapply(name,
               function(n) do.call(rbind, lapply(l, `[[`, n)) )
  names(out) = name
  return(out)
}


#### pre processing ####

# treatment variables
A = as.logical(probitsim$a2)
probitsim$a1 = NULL # unused
probitsim$a2 = NULL # unused
probitsim$a3 = NULL # unused
probitsim$a4 = NULL # unused

# outcome variable
Y = probitsim$wgt3
probitsim$wgt3 = NULL

# data transformation
probitsim$cage  = probitsim$age-mean(probitsim$age)
probitsim$cage2 = probitsim$cage^2


X = model.matrix(A ~ -1 + factor(location) + factor(educ) + cage + cage2 + factor(smoke) + factor(allergy),
                 data = probitsim)
X = matrix(data = X, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))
X = X[ ,"factor(location)4" != colnames(X)] # location has 4 levels so only needs 3 dummies

bin_var = apply(X, 2, function(x) length(unique(x))<3)
X[ ,!bin_var] = scale(X[ ,!bin_var])

# response surface
rs = lapply(list("y0"=!A, "y1"=A), function(subset)
  lm(Y ~ X, subset = subset)) |>
  lapply(predict, newdata=probitsim, type="response") |>
  lapply(as.vector)

outcome_block = outcome_blocks_linear(Y = Y, A = A, X = X)


# IPTW
iptw = IPTW_logreg_weight(X = as.data.frame(X), A = A)
ps_block_iptw = propensity_blocks_logit_mle(A = A, X = X, model.fit = iptw$model)
iptw_wls_ate  = wls_ci(Y, A, ps_block = ps_block_iptw, target = "ATE", alpha = 0.05)
iptw_wls_att  = wls_ci(Y, A, ps_block = ps_block_iptw, target = "ATT", alpha = 0.05)
iptw_dr_ate   = dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_iptw, target = "ATE", alpha = 0.05)
iptw_dr_att   = dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_iptw, target = "ATT", alpha = 0.05)


# EB
eb = EB(X = `colnames<-`(X,1:ncol(X)), A = A)
eb_ate = pmax(eb$ATE, 0); eb_att = pmax(eb$ATT, 0)
eb_wls_ate  = wls_ci_optimist(Y, A, W = eb_ate, alpha = 0.05)
eb_wls_att  = wls_ci_optimist(Y, A, W = eb_att, alpha = 0.05)
eb_dr_ate   = dr_ci_optimist(A = A, outcome_block = outcome_block, W = eb_ate, target = "ATE", alpha = 0.05)
eb_dr_att   = dr_ci_optimist(A = A, outcome_block = outcome_block, W = eb_att, target = "ATT", alpha = 0.05)


# KOM
kom.hp = KOM_hp(X = X, A = A, Y = Y)
kom = KOM_weight(X = X, A = A, lambda = lapply(kom.hp, `[[`, "lambda"), scale  = lapply(kom.hp, `[[`, "scale"))

kom_ate = pmax(kom$ATE,0); kom_att = pmax(kom$ATT, 0)
kom_wls_ate  = wls_ci_optimist(Y, A, W = kom_ate, alpha = 0.05)
kom_wls_att  = wls_ci_optimist(Y, A, W = kom_att, alpha = 0.05)
kom_dr_ate   = dr_ci_optimist(A = A, outcome_block = outcome_block, W = kom_ate, target = "ATE", alpha = 0.05)
kom_dr_att   = dr_ci_optimist(A = A, outcome_block = outcome_block, W = kom_att, target = "ATT", alpha = 0.05)


# TLF



if(inherits(try(load("TLF_param.RData")), "try-error")) {
  D1 = as.matrix(dist(X, method = "manhattan", diag = TRUE, upper = TRUE))
  sigma = 1 / median(D1)
  
  lambdas = 10^seq(-6,2,length.out=80)
  tlf.hp = kernel.hp.gridsearch(T = A, X = X, alphas = c(-1,0),
                                betas = c(-1,-1), lambdas = lambdas, sigmas = sigma,
                                nbr_fold = 5, energy = .99, intercept = TRUE, best = FALSE,
                                verbose = FALSE, seed = 123)
  tlf.lambda = lapply(list(ATE=1, ATT=2), function(i) lambdas[which.min(log(tlf.hp[ ,i]))])
  tlf.lambda$ATT = lambdas[38] # no global min: hand-picked after checking graph

  K.eigen = eigen(exp(-sigma * D1))
  save(sigma, tlf.hp, tlf.lambda, K.eigen, file="TLF_param.RData")
}
  
tlf_out = TLF_weight(X = X, A = A, sigma = sigma, lambda = tlf.lambda, K.eigen = K.eigen, intercept = TRUE)
res_ate = tlf_out$model$ATE
res_att = tlf_out$model$ATT
ps_block_tlf_ate = propensity_blocks_cbpstlf(A = A, U = res_ate$U, eta = res_ate$eta, d = res_ate$d, lambda = tlf.lambda$ATE, target="ATE")
ps_block_tlf_att = propensity_blocks_cbpstlf(A = A, U = res_att$U, eta = res_att$eta, d = res_att$d, lambda = tlf.lambda$ATT, target="ATT")
tlf_wls_ate = wls_ci(Y, A, ps_block = ps_block_tlf_ate, target = "ATE", alpha = 0.05)
tlf_wls_att = wls_ci(Y, A, ps_block = ps_block_tlf_att, target = "ATT", alpha = 0.05)
tlf_dr_ate  = dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_tlf_ate, target = "ATE", alpha = 0.05)
tlf_dr_att  = dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_tlf_att, target = "ATT", alpha = 0.05)

data = as.list(environment())
save(data, file="probitsim_output.RData")
