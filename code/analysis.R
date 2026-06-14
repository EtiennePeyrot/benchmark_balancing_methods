rm(list = ls())
path = "C:/Users/skoua/Work/project benchmark balancing/new version/sim_2909128"
load(paste0(path, "/data.RData"))


# Dimension names ----------------------------------------------------------


pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
method_names    = dimnames(data[[1]]$estimate_arr)[[1]]
effect_names    = dimnames(data[[1]]$estimate_arr)[[2]]
estimand_names  = dimnames(data[[1]]$estimate_arr)[[3]]
estimator_names = dimnames(data[[1]]$estimate_arr)[[4]]
weight_method_names = dimnames(data[[1]]$W)[[1]]
variables_names = dimnames(data[[1]]$smd)[[1]]
group_names = dimnames(data[[1]]$ess)[[1]]

true_effect = c(none = 0, small = 0.05, high = 0.10)


# Remove abnormal values --------------------------------------------------

abnormal_values = data.frame(
  n = character(),
  prop_trt = character(),
  confdg_lvl = character(),
  method = character(),
  effect = character(),
  estimand = character(),
  estimator = character(),
  rep = integer(),
  estimate = numeric(),
  stringsAsFactors = FALSE
)

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row, "confdg_lvl"]
  
  estimate_arr = data[[n, prop_trt, confdg_lvl]]$estimate_arr
  
  estimate_values = estimate_arr[, , , , "estimate", , drop = FALSE]
  
  out_of_range = !is.na(estimate_values) & (estimate_values < -1 | estimate_values > 1)
  
  if (all(!out_of_range)) next
  
  idx = which(out_of_range, arr.ind = TRUE)
  
  for (i in 1:nrow(idx)) {
    method_i    = dimnames(estimate_arr)[[1]][idx[i, 1]]
    effect_i    = dimnames(estimate_arr)[[2]][idx[i, 2]]
    estimand_i  = dimnames(estimate_arr)[[3]][idx[i, 3]]
    estimator_i = dimnames(estimate_arr)[[4]][idx[i, 4]]
    rep_i       = idx[i, 6]
    estimate_i = estimate_values[method_i, effect_i, estimand_i, estimator_i, "estimate", rep_i]
    
    abnormal_values = rbind(
      abnormal_values,
      data.frame(n = n, prop_trt = prop_trt, confdg_lvl = confdg_lvl,
                 method = method_i, effect = effect_i, estimand = estimand_i,
                 estimator = estimator_i, rep = rep_i,
                 estimate = as.numeric(estimate_i), stringsAsFactors = FALSE)
    )
    estimate_arr[method_i, effect_i, estimand_i, estimator_i, , rep_i] = NA_real_
  }
  
  data[[n, prop_trt, confdg_lvl]]$estimate_arr = estimate_arr
}

save(abnormal_values, file = paste0(path, "/abnormal_values.RData"))

# Bias --------------------------------------------------------------------


bias = array(data = NA_real_,
             dim = c(dim(data), length(method_names), length(effect_names),
                     length(estimand_names), length(estimator_names)),
             dimnames = c(dimnames(data),
                          list(method=method_names,
                               effect=effect_names,
                               estimand=estimand_names,
                               estimator=estimator_names)))

bias_mcse = bias

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  estimate_values = data[[n, prop_trt, confdg_lvl]]$estimate_arr[, , , , "estimate", ]
  estimate_error = sweep(estimate_values, 2, true_effect, "-")

  bias[n, prop_trt, confdg_lvl, , , , ] =
    apply(estimate_error, 1:4, mean, na.rm = TRUE)

  bias_mcse[n, prop_trt, confdg_lvl, , , , ] =
    apply(estimate_error, 1:4, function(x) {
      sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
    })
}


# Variance ----------------------------------------------------------------


variance = array(data = NA_real_, dim = c(dim(data), length(method_names), length(effect_names), length(estimand_names), length(estimator_names)),
                 dimnames = c(dimnames(data),
                              list(method=method_names,
                                   effect=effect_names,
                                   estimand=estimand_names,
                                   estimator=estimator_names)))

variance_mcse = variance

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  estimate_values = data[[n, prop_trt, confdg_lvl]]$estimate_arr[, , , , "estimate", ]

  variance[n, prop_trt, confdg_lvl, , , , ] =
    apply(estimate_values, 1:4, var, na.rm = TRUE)

  variance_mcse[n, prop_trt, confdg_lvl, , , , ] =
  apply(estimate_values, 1:4, function(x) {
    x = x[!is.na(x)]
    R = length(x)
    Q = (x - mean(x))^2
    sd(Q) / sqrt(R)
  })
}


# Mean Absolute Error (MAE) -----------------------------------------------


mae = array(data = NA_real_, dim = c(dim(data), length(method_names), length(effect_names), length(estimand_names), length(estimator_names)),
            dimnames = c(dimnames(data),
                         list(method=method_names,
                              effect=effect_names,
                              estimand=estimand_names,
                              estimator=estimator_names)))

mae_mcse = mae

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  estimate_values = data[[n, prop_trt, confdg_lvl]]$estimate_arr[, , , , "estimate", ]
  estimate_error = sweep(estimate_values, 2, true_effect, "-")

  mae[n, prop_trt, confdg_lvl, , , , ] =
    apply(estimate_error, 1:4, function(x) mean(abs(x), na.rm = TRUE))

  mae_mcse[n, prop_trt, confdg_lvl, , , , ] =
    apply(estimate_error, 1:4, function(x) {
      sd(abs(x), na.rm = TRUE) / sqrt(sum(!is.na(x)))
    })
}


# Root Mean Squared Error (RMSE) ------------------------------------------


rmse = array(data = NA_real_, dim = c(dim(data), length(method_names), length(effect_names), length(estimand_names), length(estimator_names)),
             dimnames = c(dimnames(data),
                          list(method=method_names,
                               effect=effect_names,
                               estimand=estimand_names,
                               estimator=estimator_names)))

rmse_mcse = rmse

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  estimate_values = data[[n, prop_trt, confdg_lvl]]$estimate_arr[, , , , "estimate", ]
  estimate_error = sweep(estimate_values, 2, true_effect, "-")

  rmse[n, prop_trt, confdg_lvl, , , , ] =
    apply(estimate_error, 1:4, function(x) sqrt(mean(x^2, na.rm = TRUE)))

  rmse_mcse[n, prop_trt, confdg_lvl, , , , ] =
    apply(estimate_error, 1:4, function(x) {
      m = mean(x^2, na.rm = TRUE)
      sd(x^2, na.rm = TRUE) / sqrt(sum(!is.na(x))) / (2 * sqrt(m))
    })
}


# CI coverage -------------------------------------------------------------


CI_coverage = array(data = NA_real_, dim = c(dim(data), length(method_names), length(effect_names), length(estimand_names), length(estimator_names)),
                    dimnames = c(dimnames(data),
                                 list(method=method_names,
                                      effect=effect_names,
                                      estimand=estimand_names,
                                      estimator=estimator_names)))

CI_coverage_mcse = CI_coverage

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  for (effect in effect_names) {
    CI_coverage[n, prop_trt, confdg_lvl, , effect, , ] =
      apply(data[[n, prop_trt, confdg_lvl]]$estimate_arr[, effect, , , , ], 1:3, function(x) {
        out = (true_effect[effect] >= x["ci_lower",]) &  (true_effect[effect] <= x["ci_upper",])
        mean(out, na.rm = TRUE)
      })

    CI_coverage_mcse[n, prop_trt, confdg_lvl, , effect, , ] =
      apply(data[[n, prop_trt, confdg_lvl]]$estimate_arr[, effect, , , , ], 1:3, function(x) {
        out = (true_effect[effect] >= x["ci_lower",]) &  (true_effect[effect] <= x["ci_upper",])
        p = mean(out, na.rm = TRUE)
        sqrt(p * (1 - p) / sum(!is.na(out)))
      })
  }
}


# Variance ratio ----------------------------------------------------------


variance_ratio = array(data = NA_real_, dim = c(dim(data), length(method_names), length(effect_names), length(estimand_names), length(estimator_names)),
                       dimnames = c(dimnames(data),
                                    list(method=method_names,
                                         effect=effect_names,
                                         estimand=estimand_names,
                                         estimator=estimator_names)))

variance_ratio_mcse = variance_ratio

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  variance_ratio[n, prop_trt, confdg_lvl, , , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$estimate_arr,1:4, function(x) {
      estimate = x["estimate",]
      estimated_variance = x["se",]^2
      ok = !is.na(estimate) & !is.na(estimated_variance)
      mean(estimated_variance[ok]) / var(estimate[ok])
    })

  variance_ratio_mcse[n, prop_trt, confdg_lvl, , , , ] =
  apply(data[[n, prop_trt, confdg_lvl]]$estimate_arr,1:4, function(x) {
    estimate = x["estimate",]
    estimated_variance = x["se",]^2
    ok = !is.na(estimate) & !is.na(estimated_variance)

    estimate = estimate[ok]
    estimated_variance = estimated_variance[ok]

    A = mean(estimated_variance)
    B = var(estimate)

    Q = (estimate - mean(estimate))^2
    Q_bar = mean(Q)

    phi_variance_ratio = (estimated_variance - A) / B -
      A * (Q - Q_bar) / B^2

    sd(phi_variance_ratio) / sqrt(length(estimate))
  })
}


# Missing values ----------------------------------------------------------


missing_value = function(name) {
  out = array(data = NA_real_, dim = c(dim(data), length(method_names), length(effect_names), length(estimand_names), length(estimator_names)),
              dimnames = c(dimnames(data),
                           list(method=method_names,
                                effect=effect_names,
                                estimand=estimand_names,
                                estimator=estimator_names)))

  for (row in 1:nrow(pars)) {
    n          = pars[row, "n"]
    prop_trt   = pars[row, "prop_trt"]
    confdg_lvl = pars[row , "confdg_lvl"]

    out[n, prop_trt, confdg_lvl, , , , ] =
      apply(data[[n, prop_trt, confdg_lvl]]$estimate_arr[, , , , name, ],1:4, function(x) mean(is.na(x)) )
  }

  return(out)
}

missing = list(
  estimate = missing_value("estimate"),
  standard_error = missing_value("se")
)

metric_missingness = list(
  bias = "estimate",
  variance = "estimate",
  mae = "estimate",
  rmse = "estimate",
  CI_coverage = c("estimate", "standard_error"),
  variance_ratio = c("estimate", "standard_error")
)


# min sum W<0 -------------------------------------------------------------


SNegW = array(data = NA_real_, dim = c(dim(data), length(weight_method_names), length(estimand_names)),
              dimnames = c(dimnames(data),
                           list(method=weight_method_names,
                                estimand=estimand_names)))

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  SNegW[n, prop_trt, confdg_lvl, , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$W,c(1,2,4), function(w) {
      sum(pmin(w,0))
    }) |> apply(1:2, min)
}


# Standardized Mean Difference (SMD) --------------------------------------


smd = array(data = NA_real_, dim = c(dim(data), length(variables_names), length(weight_method_names), length(estimand_names)),
            dimnames = c(dimnames(data),
                         list(variables=variables_names,
                              method=weight_method_names,
                              estimand=estimand_names)))

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  smd[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$smd,1:3, mean, na.rm = TRUE)
}


# Effective Sample Size (ESS) ---------------------------------------------


ess = array(data = NA_real_, dim = c(dim(data), length(group_names), length(weight_method_names), length(estimand_names)),
            dimnames = c(dimnames(data),
                         list(group=group_names,
                              method=weight_method_names,
                              estimand=estimand_names)))

for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]

  ess[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$ess,1:3, median, na.rm = TRUE)
}


save(bias, bias_mcse,
     variance, variance_mcse,
     mae, mae_mcse,
     rmse, rmse_mcse,
     CI_coverage, CI_coverage_mcse,
     variance_ratio, variance_ratio_mcse,
     missing, metric_missingness,
     SNegW, smd, ess,
     file = paste0(path,"/metrics.RData"))
