rm(list = ls())
graphics.off()
set.seed(123)

# import
source("TLF_weight.R")
source("../GenerateData.R")


task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ntasks  = as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
outdir = as.character(Sys.getenv("OUTDIR"))
N_rep = as.integer(Sys.getenv("N_REPEAT"))


idx = which(rep(1:ntasks, length.out = N_rep) == task_id)
N = length(idx) # number of iteration for current task

sigmas_prod = 2^(-15:1) # search around sigma heuristic
lambdas = c(0, 10^(-7:2))

# function

norm = function(w,a) w / ifelse(a, sum(w[a]), sum(w[!a]))
smd_weighted = function(X, A, W, sddenom = c("pooled","treated","control","all")) {
  sddenom = match.arg(sddenom)
  W0 = W * !A; W1 = W * A
  m0 = sum(X * W0) / sum(W0); m1 = sum(X * W1) / sum(W1)
  v0 = sum(W0 * (X - m0)^2) / sum(W0); v1 = sum(W1 * (X - m1)^2) / sum(W1)
  
  sd_den = switch(sddenom,
                  pooled  = sqrt((v1 + v0) / 2),
                  treated = sqrt(v1),
                  control = sqrt(v0),
                  all     = sqrt(sum(W * (X - sum(W * X) / sum(W))^2) / sum(W))
  )
  return((m1 - m0) / sd_den)
}

max_smd = function(X, A, W, sddenom) {
  max(abs(sapply(X,smd_weighted,A=A, W=W, sddenom=sddenom)))
}

cv = function(W) sd(W) / sum(W)


# core
res = array(data = NA_real_,
            dim = c(2,4,3,3,length(lambdas),length(sigmas_prod),5,N),
            dimnames = list(estimand=c("ATE", "ATT"),
                            n=c(250, 500, 1000, 2000),
                            prop_trt=c("low", "moderate", "high"),
                            confdg_lvl=c("low", "moderate", "high"),
                            lambda=lambdas,
                            sigma=paste0("sig0 2^", log2(sigmas_prod)),
                            metric=c("grad CBSR", "max SMD", "weight CV", "Effect", "Effect norm"),
                            seed = 1:N))


pars1 = do.call(expand.grid,c(dimnames(res)[c(2:4,8)], stringsAsFactors = FALSE))
pars2 = expand.grid(lambda_idx=seq_along(lambdas), sigma_idx=seq_along(sigmas_prod))

for(row1 in 1:nrow(pars1)) {
  n = pars1[row1,"n"]; prop_trt = pars1[row1,"prop_trt"]
  confdg_lvl = pars1[row1, "confdg_lvl"]; i = as.integer(pars1[row1,"seed"])
  
  set.seed(idx[i])
  data = GenerateData(n=as.integer(n), prop_trt=prop_trt, confdg_lvl=confdg_lvl)
  
  X = data$X
  A = as.logical(data$A)
  Y = data$Y
  D = as.matrix(dist(X, "manhattan"))
  
  sigma_heur = 1 / median(D)
  sigmas = sigma_heur * sigmas_prod
  
  try({
    res["ATE", n, prop_trt, confdg_lvl, , , "grad CBSR", i] =
      kernel.hp.gridsearch(T = A, X = X, alphas = -1, betas = -1,
                           lambdas = lambdas, sigmas = sigmas,
                           energy = .9, intercept = TRUE,
                           nbr_fold = 5, best = F)
  })
  try({
    res["ATT", n, prop_trt, confdg_lvl, , , "grad CBSR", i] =
      kernel.hp.gridsearch(T = A, X = X, alphas = 0, betas = -1,
                           lambdas = lambdas, sigmas = sigmas,
                           energy = .9, intercept = TRUE,
                           nbr_fold = 5, best = F)
  })
  
  sigma_ = max(sigmas) + 1
  for (row2 in 1:nrow(pars2)) try({
    lambda_idx = pars2[row2, 1]; lambda = lambdas[lambda_idx]
    sigma_idx = pars2[row2, 2]; sigma = sigmas[sigma_idx]
    
    # compute eigen decomposition iif sigma changes
    if (sigma_ != sigma) {
      sigma_ = sigma
      K = exp(-sigma * D)
      K.eigen = eigen(K)
    }
    
    W = TLF_weight(X = X, A = A, lambda = list(ATE=lambda, ATT=lambda),
                   K.eigen = K.eigen, intercept = TRUE,
                   maxiter = 100, tol = 1e-4, energy = .9,
                   verbose = FALSE)$weights
    
    res["ATE", n, prop_trt, confdg_lvl, lambda_idx, sigma_idx, "max SMD", i] =
      max_smd(X = X, A = A, W = W$ATE, sddenom = "all")
    res["ATT", n, prop_trt, confdg_lvl, lambda_idx, sigma_idx, "max SMD", i] =
      max_smd(X = X, A = A, W = W$ATT, sddenom = "all")
    
    res["ATE", n, prop_trt, confdg_lvl, lambda_idx, sigma_idx, "weight CV", i] =
      cv(W = W$ATE)
    res["ATT", n, prop_trt, confdg_lvl, lambda_idx, sigma_idx, "weight CV", i] =
      cv(W = W$ATT)
    
    res["ATE", n, prop_trt, confdg_lvl, lambda_idx, sigma_idx, "Effect", i] =
      sum(W$ATE * (2 * A - 1) * Y)
    res["ATT", n, prop_trt, confdg_lvl, lambda_idx, sigma_idx, "Effect", i] =
      sum(W$ATT * (2 * A - 1) * Y)
    
    res["ATE", n, prop_trt, confdg_lvl, lambda_idx, sigma_idx, "Effect norm", i] =
      sum(norm(W$ATE, A) * (2 * A - 1) * Y)
    res["ATT", n, prop_trt, confdg_lvl, lambda_idx, sigma_idx, "Effect norm", i] =
      sum(norm(W$ATT, A) * (2 * A - 1) * Y)
  })
}


save(res, file=paste0(outdir,"/res/",task_id,".RData"))