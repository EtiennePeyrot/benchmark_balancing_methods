rm(list = ls())

resdir = as.character(Sys.getenv("RESDIR"))

n_values = c(250, 500, 1000, 2000)
prop_trt_values = c("low", "moderate", "high")
confdg_lvl_values = c("low", "moderate", "high")
NREP = 5000

methods = c("IPTW", "CBPS-JI", "CBPS-OI", "CBPS-TLF", "EB", "KOM")
effects = c("none", "small", "high")
weight_methods = c("IPTW", "CBPS-JI", "CBPS-OI", "CBPS-TLF", "EB", "KOM_none", "KOM_small", "KOM_high")
estimands = c("ATE", "ATT")
estimators = c("WLS", "DR")
values = c("estimate", "ci_lower", "ci_upper", "se")
groups = c("treated", "control")
variables = paste0("V", 1:10)
seeds = as.character(1:NREP)

out = array(data = list(list()),
            dim = c(length(n_values), length(prop_trt_values), length(confdg_lvl_values)),
            dimnames = list(n = n_values,
                            prop_trt = prop_trt_values,
                            confdg_lvl = confdg_lvl_values))

# X stored only by n, since X depends only on n and seed
X = setNames(vector("list", length(n_values)), as.character(n_values))

pars = expand.grid(seed = 1:NREP,
                   n = n_values,
                   prop_trt = prop_trt_values,
                   confdg_lvl = confdg_lvl_values,
                   stringsAsFactors = FALSE)

for (row in 1:nrow(pars)) {
  # cat("\r", format(round(100 * row / nrow(pars), 1), width = 5, justify = "right"), "%   ")
  
  n_val = pars[row, "n"]
  prop_trt = pars[row, "prop_trt"]
  confdg_lvl = pars[row, "confdg_lvl"]
  seed = pars[row, "seed"]
  
  file_path = paste0(resdir, "/", n_val, "obs_", prop_trt, "proptrt_",
                     confdg_lvl, "confdglvl_", seed, ".RData")
  is_loaded = try(load(file_path))
  if (inherits(is_loaded, "try-error")) next
  
  Y_mat = cbind(none = data$Y$none,
                small = data$Y$small,
                high = data$Y$high)
  
  # init per-setting container, without X
  if (seed == 1) {
    out[[as.character(n_val), prop_trt, confdg_lvl]] = list(
      A = array(data = NA_real_,
                dim = c(n_val, NREP),
                dimnames = list(unit = NULL, seed = seeds)),
      W = array(data = NA_real_,
                dim = c(length(weight_methods), length(estimands), n_val, NREP),
                dimnames = list(method = weight_methods,
                                estimand = estimands,
                                unit = NULL,
                                seed = seeds)),
      Y = array(data = NA_real_,
                dim = c(n_val, length(effects), NREP),
                dimnames = list(unit = NULL,
                                effect = effects,
                                seed = seeds)),
      smd = array(data = NA_real_,
                  dim = c(length(variables), length(weight_methods), length(estimands), NREP),
                  dimnames = list(variable = variables,
                                  method = weight_methods,
                                  estimand = estimands,
                                  seed = seeds)),
      ess = array(data = NA_real_,
                  dim = c(length(groups), length(weight_methods), length(estimands), NREP),
                  dimnames = list(group = groups,
                                  method = weight_methods,
                                  estimand = estimands,
                                  seed = seeds)),
      estimate_arr = array(data = NA_real_,
                           dim = c(length(methods), length(effects), length(estimands), length(estimators), length(values), NREP),
                           dimnames = list(method = methods,
                                           effect = effects,
                                           estimand = estimands,
                                           estimator = estimators,
                                           value = values,
                                           seed = seeds))
    )
  }
  
  # init per-n X container on first seed
  if (seed == 1 && length(X[[as.character(n_val)]]) == 0) {
    X[[as.character(n_val)]] = array(data = NA_real_,
                                     dim = c(n_val, length(variables), NREP),
                                     dimnames = list(unit = NULL,
                                                     variable = variables,
                                                     seed = seeds))
  }
  
  # fill per-setting objects
  out[[as.character(n_val), prop_trt, confdg_lvl]]$A[ , seed] = data$A
  out[[as.character(n_val), prop_trt, confdg_lvl]]$W[ , , , seed] = W
  out[[as.character(n_val), prop_trt, confdg_lvl]]$Y[ , , seed] = Y_mat
  out[[as.character(n_val), prop_trt, confdg_lvl]]$smd[ , , , seed] = smd
  out[[as.character(n_val), prop_trt, confdg_lvl]]$ess[ , , , seed] = ess
  out[[as.character(n_val), prop_trt, confdg_lvl]]$estimate_arr[ , , , , , seed] = estimate_arr
  
  # fill per-n X
  X[[as.character(n_val)]][ , , seed] = as.matrix(data$X)
}

data = out
rm("out")
save(data, X, file = "data.RData")
