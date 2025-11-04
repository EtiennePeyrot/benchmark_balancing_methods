rm(list = ls())
path = "C:/Users/skoua/Work/project benchmark balancing/new version/sim_2284941/"
load(paste0(path,"data.RData"))

# Bias --------------------------------------------------------------------


bias = array(data = NA_real_, dim = c(dim(data),4,2,2),
             dimnames = c(dimnames(data),
                          list(method=c("IPTW", "EB", "KOM", "TLF"),
                               estimand=c("ATE", "ATT"),
                               estimator=c("WLS", "DR"))))
pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  bias[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$estimate[,,,"estimate",],1:3, mean, na.rm = T)
}


# Variance ----------------------------------------------------------------


variance = array(data = NA_real_, dim = c(dim(data),4,2,2),
                 dimnames = c(dimnames(data),
                              list(method=c("IPTW", "EB", "KOM", "TLF"),
                                   estimand=c("ATE", "ATT"),
                                   estimator=c("WLS", "DR"))))
pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  variance[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$estimate[,,,"estimate",],1:3, var, na.rm = T)
}


# Mean Absolute Error (MAE) -----------------------------------------------


mae = array(data = NA_real_, dim = c(dim(data),4,2,2),
                 dimnames = c(dimnames(data),
                              list(method=c("IPTW", "EB", "KOM", "TLF"),
                                   estimand=c("ATE", "ATT"),
                                   estimator=c("WLS", "DR"))))
pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  mae[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$estimate[,,,"estimate",],1:3, function(x) mean(abs(x), na.rm=T))
}



# Residual Mean Squared Error (RMSE) --------------------------------------


rmse = array(data = NA_real_, dim = c(dim(data),4,2,2),
            dimnames = c(dimnames(data),
                         list(method=c("IPTW", "EB", "KOM", "TLF"),
                              estimand=c("ATE", "ATT"),
                              estimator=c("WLS", "DR"))))
pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  rmse[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$estimate[,,,"estimate",],1:3, function(x) mean(x^2, na.rm=T))
}


# CI coverage -------------------------------------------------------------


CI_coverage = array(data = NA_real_, dim = c(dim(data),4,2,2),
                    dimnames = c(dimnames(data),
                                 list(method=c("IPTW", "EB", "KOM", "TLF"),
                                      estimand=c("ATE", "ATT"),
                                      estimator=c("WLS", "DR"))))
pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  CI_coverage[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$estimate,1:3, function(x) {
      out = (0 >= x["ci_lower",]) &  (0 <= x["ci_upper",])
      sapply(out, isTRUE) |> mean()
    })
}


# CI width ----------------------------------------------------------------


CI_width = array(data = NA_real_, dim = c(dim(data),4,2,2),
                    dimnames = c(dimnames(data),
                                 list(method=c("IPTW", "EB", "KOM", "TLF"),
                                      estimand=c("ATE", "ATT"),
                                      estimator=c("WLS", "DR"))))
pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  CI_width[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$estimate,1:3, function(x) {
      mean((x["ci_upper",] - x["ci_lower",]), na.rm=T)
    })
}



# Failure rate ------------------------------------------------------------


fail = array(data = NA_real_, dim = c(dim(data),4,2,2),
                 dimnames = c(dimnames(data),
                              list(method=c("IPTW", "EB", "KOM", "TLF"),
                                   estimand=c("ATE", "ATT"),
                                   estimator=c("WLS", "DR"))))
pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  fail[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$estimate[,,,"estimate",],1:3, function(x) {
      mean(is.na(x))
    })
}


# min sum W<0 -------------------------------------------------------------


SNegW = array(data = NA_real_, dim = c(dim(data),4,2),
             dimnames = c(dimnames(data),
                          list(method=c("IPTW", "EB", "KOM", "TLF"),
                               estimand=c("ATE", "ATT"))))
pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  SNegW[n, prop_trt, confdg_lvl, , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$W,c(1,2,4), function(w) {
      sum(pmin(w,0))
    }) |> apply(1:2, min)
}


# Standardize Mean Difference (SMD) ---------------------------------------


smd = array(data = NA_real_, dim = c(dim(data),10,4,2),
            dimnames = c(dimnames(data),
                         list(variables=paste0("X",1:10),
                              method=c("IPTW", "EB", "KOM", "TLF"),
                              estimand=c("ATE", "ATT"))))

pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  smd[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$smd,1:3, mean, na.rm = T)
}


# Effective Sample Size (ESS) ---------------------------------------------


ess = array(data = NA_real_, dim = c(dim(data),2,4,2),
            dimnames = c(dimnames(data),
                         list(group=c("treated", "control"),
                              method=c("IPTW", "EB", "KOM", "TLF"),
                              estimand=c("ATE", "ATT"))))

pars = do.call(expand.grid, c(dimnames(data), stringsAsFactors=F))
for (row in 1:nrow(pars)) {
  n          = pars[row, "n"]
  prop_trt   = pars[row, "prop_trt"]
  confdg_lvl = pars[row , "confdg_lvl"]
  
  ess[n, prop_trt, confdg_lvl, , , ] =
    apply(data[[n, prop_trt, confdg_lvl]]$ess,1:3, median, na.rm = T)
}



save(bias, variance, mae, rmse, CI_coverage, CI_width, fail, SNegW, smd, ess,
     file = paste0(path,"metrics.RData"))
