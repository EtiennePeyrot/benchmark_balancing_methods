rm(list = ls())

resdir = as.character(Sys.getenv("RESDIR"))

n = c(250,500,1000,2000)
prop_trt = c("low", "moderate", "high")
confdg_lvl = c("low", "moderate", "high")
NREP = 5000

out = array(data = list(list()),
            dim = sapply(list(n, prop_trt, confdg_lvl), length),
            dimnames = list(n=n, prop_trt=prop_trt, confdg_lvl=confdg_lvl))

# X stored only by n (since X depends only on n and seed)
X = setNames(vector("list", length(n)), as.character(n))

pars = do.call(expand.grid, c(list(seed=1:NREP),dimnames(out), stringsAsFactors=F))

for (row in 1:nrow(pars)) {
  # cat("\r", format(round(row / 1800,1), width=3, justify="right"), "%   ")
  
  n = pars[row, "n"]; confdg_lvl = pars[row, "confdg_lvl"]
  prop_trt = pars[row, "prop_trt"]; seed = pars[row, "seed"]
  file_path = paste0(resdir,"/",n,"obs_",prop_trt,"proptrt_",confdg_lvl,"confdglvl_",seed,".RData")
  is_loaded = try(load(file_path))
  if (inherits(is_loaded, "try-error")) next
  
  # init per-setting container (without X)
  if (seed==1) out[[n,prop_trt,confdg_lvl]] = list(
    A = array(data = NA_real_, dim = c(length(data$A), NREP)),
    W = array(data = NA_real_, dim = c(dim(W), NREP), dimnames = dimnames(W)),
    Y = array(data = NA_real_, dim = c(length(data$Y), NREP)),
    smd = array(data = NA_real_, dim = c(dim(smd), NREP), dimnames = dimnames(smd)),
    ess = array(data = NA_real_, dim = c(dim(ess), NREP), dimnames = dimnames(ess)),
    estimate_arr = array(data = NA_real_,
                         dim = c(dim(estimate_arr)[1:4], NREP),
                         dimnames = dimnames(estimate_arr)[1:4])
  )
  
  # init per-n X container on first encounter
  if (seed==1 && length(X[[as.character(n)]])==0) {
    X[[as.character(n)]] = array(data = NA_real_,
                                 dim = c(dim(as.matrix(data$X)), NREP),
                                 dimnames = dimnames(data$X))
  }
  
  # fill per-setting objects
  out[[n,prop_trt,confdg_lvl]]$A[ ,seed] = data$A
  out[[n,prop_trt,confdg_lvl]]$W[ , , ,seed] = W
  out[[n,prop_trt,confdg_lvl]]$Y[ ,seed] = data$Y
  out[[n,prop_trt,confdg_lvl]]$smd[ , , ,seed] = smd
  out[[n,prop_trt,confdg_lvl]]$ess[ , , ,seed] = ess
  out[[n,prop_trt,confdg_lvl]]$estimate_arr[ , , , ,seed] = estimate_arr # for old version: estimate_arr[ , , , ,1,drop=F]
  
  # fill per-n X
  X[[as.character(n)]][ , , seed] = as.matrix(data$X)
}

data = out; rm("out")
save(data, X, file = "data.RData")
