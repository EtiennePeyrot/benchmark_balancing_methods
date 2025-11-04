rm(list = ls())
graphics.off()
set.seed(123)

# import
source("GenerateData.R")
source("IPTW/IPTW_weight.R")
source("IPTW/confidence_interval.R")
source("EB/EB_weight.R")
source("KOM/KOM_hp_selec.R")
source("KOM/KOM_weight.R")
source("KOM/confidence_interval.R")
source("TLF/TLF_weight.R")
source("TLF/confidence_interval.R")


task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ntasks  = as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
outdir = as.character(Sys.getenv("OUTDIR"))
N_rep = as.integer(Sys.getenv("N_REPEAT"))


idx = which(rep(1:ntasks, length.out = N_rep) == task_id)
N = length(idx) # number of iteration for current task

# helper functions
smd_weighted = function(X, A, w) {
  A  = as.integer(A)
  wt = w * A
  wc = w * (1 - A)
  X  = as.matrix(X)
  mt = colSums(X * wt) / sum(wt)
  mc = colSums(X * wc) / sum(wc)
  Xt = sweep(X, 2, mt, "-"); Xc = sweep(X, 2, mc, "-")
  vt = colSums((Xt^2) * (wt / sum(wt)))
  vc = colSums((Xc^2) * (wc / sum(wc)))
  (mt - mc) / sqrt((vt + vc) / 2)
}

ess_by_group = function(w, A) {
  A  = as.integer(A)
  wt = w * A; wc = w * (1 - A)
  c(treated = (sum(wt)^2) / sum(wt^2),
    control = (sum(wc)^2) / sum(wc^2))
}

unpack = function(wls_ate, wls_att, dr_ate, dr_att) {
  estimands  = c("ATE","ATT")
  estimators = c("WLS","AIPW")
  values     = c("estimate","ci_lower","ci_upper","se")
  out = array(NA_real_, dim=c(length(estimands), length(estimators), length(values)),
              dimnames=list(estimand=estimands, estimator=estimators, value=values))
  try(out["ATE","WLS",]  <- c(wls_ate$tau, wls_ate$ci[1], wls_ate$ci[2], wls_ate$se))
  try(out["ATT","WLS",]  <- c(wls_att$tau, wls_att$ci[1], wls_att$ci[2], wls_att$se))
  try(out["ATE","AIPW",] <- c(dr_ate$tau,  dr_ate$ci[1],  dr_ate$ci[2],  dr_ate$se))
  try(out["ATT","AIPW",] <- c(dr_att$tau,  dr_att$ci[1],  dr_att$ci[2],  dr_att$se))
  out
}

# pre computed hyper parameters
tlf.lambda = array(
  data = c(0.01 , 0.001, 0.01 , 0.001, 0.001, 1e-04, 0.001, 1e-04, 0.01 ,
           0.01 , 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.01 , 0.01 ,
           0.01 , 0.01 , 0.01 , 0.01 , 0.001, 0.001, 0.01 , 0.001, 0.01 ,
           0.001, 0.01 , 0.001, 0.01 , 0.001, 0.01 , 0.01 , 0.01 , 0.01 ,
           0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 ,
           0.01 , 0.01 , 0.01 , 0.1  , 0.01 , 0.01 , 0.001, 0.01 , 0.001,
           0.01 , 0.001, 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 ,
           0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 , 0.01 ),
  dim = c(2,4,3,3),
  dimnames = list(estimand=c("ATE", "ATT"),
                  n=c(250, 500, 1000, 2000),
                  prop_trt=c("low", "moderate", "high"),
                  confdg_lvl=c("low", "moderate", "high")))



# core

pars = expand.grid(prop_trt = c("low", "moderate", "high"),
                   confdg_lvl = c("low", "moderate", "high"),
                   n = c(250, 500, 1000, 2000),
                   seed_idx = 1:N, stringsAsFactors = F)

X_ = 0
K.eigen_ = 0
for(row in 1:nrow(pars)) {
  n = pars[row, "n"]; i = as.integer(pars[row, "seed_idx"])
  confdg_lvl = pars[row, "confdg_lvl"]; prop_trt = pars[row, "prop_trt"]
  
  set.seed(idx[i])
  
  data = GenerateData(n=as.integer(n), prop_trt=prop_trt, confdg_lvl=confdg_lvl)
  X = data$X
  A = as.logical(data$A)
  Y = data$Y
  
  methods   = c("IPTW.logreg","EB","KOM","TLF")
  estimands = c("ATE","ATT")
  estimators = c("WLS","AIPW")
  values     = c("estimate","ci_lower","ci_upper","se")
  seed_lab   = as.character(idx[i])
  
  W = array(NA_real_, dim=c(length(methods), length(estimands), as.integer(n)),
            dimnames=list(method=methods, estimand=estimands))
  
  smd = array(NA_real_, dim=c(ncol(X), length(methods), length(estimands)),
              dimnames=list(variable=colnames(X), method=methods, estimand=estimands))
  ess = array(NA_real_, dim=c(2, length(methods), length(estimands)),
              dimnames=list(group=c("treated","control"), method=methods, estimand=estimands))
  
  estimate_arr = array(NA_real_, dim=c(length(methods), length(estimands), length(estimators), length(values)),
                       dimnames=list(method=methods, estimand=estimands, estimator=estimators, value=values))
  
  outcome_block = outcome_blocks_binomial(Y = Y, A = A, X = X)
  
  
  # IPTW
  iptw = try(IPTW_logreg_weight(X = X, A = A))
  try(W["IPTW.logreg", , ] <- do.call(rbind, iptw$W))
  ps_block_iptw = try(propensity_blocks_logit_mle(A = A, X = X, model.fit = iptw$model))
  iptw_wls_ate  = try(wls_ci(Y, A, ps_block = ps_block_iptw, target = "ATE", alpha = 0.05))
  iptw_wls_att  = try(wls_ci(Y, A, ps_block = ps_block_iptw, target = "ATT", alpha = 0.05))
  iptw_dr_ate   = try(dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_iptw, target = "ATE", alpha = 0.05))
  iptw_dr_att   = try(dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_iptw, target = "ATT", alpha = 0.05))
  estimate_arr["IPTW.logreg", , , ] = unpack(iptw_wls_ate, iptw_wls_att, iptw_dr_ate, iptw_dr_att)
  
  
  # EB
  eb = try(EB(X = X, A = A))
  try(W["EB", , ] <- do.call(rbind, eb))
  eb_ate = try(pmax(eb$ATE, 0)); eb_att = try(pmax(eb$ATT, 0))
  eb_wls_ate  = try(wls_ci_optimist(Y, A, W = eb_ate, alpha = 0.05))
  eb_wls_att  = try(wls_ci_optimist(Y, A, W = eb_att, alpha = 0.05))
  eb_dr_ate   = try(dr_ci_optimist(A = A, outcome_block = outcome_block, W = eb_ate, target = "ATE", alpha = 0.05))
  eb_dr_att   = try(dr_ci_optimist(A = A, outcome_block = outcome_block, W = eb_att, target = "ATT", alpha = 0.05))
  estimate_arr["EB", , , ] = unpack(eb_wls_ate, eb_wls_att, eb_dr_ate, eb_dr_att)
  
  
  # KOM
  kom.hp = try(KOM_hp(X = X, A = A, Y = Y))
  kom = try(KOM_weight(X = X, A = A, lambda = lapply(kom.hp, `[[`, "lambda"), scale  = lapply(kom.hp, `[[`, "scale")))
  try(W["KOM", , ] <- do.call(rbind, kom))
  kom_ate = try(pmax(kom$ATE,0)); kom_att <- try(pmax(kom$ATT, 0))
  kom_wls_ate  = try(wls_ci_optimist(Y, A, W = kom_ate, alpha = 0.05))
  kom_wls_att  = try(wls_ci_optimist(Y, A, W = kom_att, alpha = 0.05))
  kom_dr_ate   = try(dr_ci_optimist(A = A, outcome_block = outcome_block, W = kom_ate, target = "ATE", alpha = 0.05))
  kom_dr_att   = try(dr_ci_optimist(A = A, outcome_block = outcome_block, W = kom_att, target = "ATT", alpha = 0.05))
  estimate_arr["KOM", , , ] = unpack(kom_wls_ate, kom_wls_att, kom_dr_ate, kom_dr_att)
  
  
  # TLF
  if (nrow(X) != nrow(X_) || !all(X == X_)) {
    D1 = as.matrix(dist(X, method = "manhattan", diag = TRUE, upper = TRUE))
    sigma = 1 / median(D1)
    K.eigen = eigen(exp(-sigma * D1))
    X_ = X; K.eigen_ = K.eigen
  } else {
    K.eigen = K.eigen_
  }
  lambda_vec = try(as.list(tlf.lambda[ ,as.character(n), prop_trt, confdg_lvl]))
  tlf_out = try(TLF_weight(X = X, A = A, sigma = sigma, lambda = lambda_vec, K.eigen = K.eigen, intercept = TRUE))
  try(W["TLF", , ] <- do.call(rbind, tlf_out$weights))
  res_ate = try(tlf_out$model$ATE)
  res_att = try(tlf_out$model$ATT)
  ps_block_tlf_ate = try(propensity_blocks_cbpstlf(A = A, U = res_ate$U, eta = res_ate$eta, d = res_ate$d, lambda = lambda_vec$ATE, target="ATE"))
  ps_block_tlf_att = try(propensity_blocks_cbpstlf(A = A, U = res_att$U, eta = res_att$eta, d = res_att$d, lambda = lambda_vec$ATT, target="ATT"))
  tlf_wls_ate = try(wls_ci(Y, A, ps_block = ps_block_tlf_ate, target = "ATE", alpha = 0.05))
  tlf_wls_att = try(wls_ci(Y, A, ps_block = ps_block_tlf_att, target = "ATT", alpha = 0.05))
  tlf_dr_ate  = try(dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_tlf_ate, target = "ATE", alpha = 0.05))
  tlf_dr_att  = try(dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_tlf_att, target = "ATT", alpha = 0.05))
  estimate_arr["TLF", , , ] = unpack(tlf_wls_ate, tlf_wls_att, tlf_dr_ate, tlf_dr_att)
  
  # SMD & ESS from W
  for (m in methods) for (tgt in estimands) {
    w = pmax(as.numeric(W[m, tgt, ]), 0)
    try(smd[, m, tgt] <- smd_weighted(X, A, w))
    try(ess[, m, tgt] <- ess_by_group(w, A))
  }
  
  save(data, W, smd, ess, estimate_arr,
       file=paste0(outdir,"/res/",n,"obs_",prop_trt,"proptrt_",confdg_lvl,"confdglvl_",idx[i],".RData"))
}
