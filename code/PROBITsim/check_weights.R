#### initialization ####

setwd("C:/Users/skoua/project benchmark balancing/PROBITsim")
rm(list = ls())
graphics.off()
set.seed(123)

# library import

library(foreign) # import dta file
library(randomForest) # random forest
library(caret) # trainControl
# library(caret)

# data import

probitsim = as.data.frame(read.dta("PROBITsim2018_v12.dta", convert.factors = FALSE))


#### pre processing ####

# treatment variables
a1 = as.logical(probitsim$a1); probitsim$a1 = NULL
a2 = as.logical(probitsim$a2); probitsim$a2 = NULL
a3 = as.logical(probitsim$a3); probitsim$a3 = NULL
a4 = as.logical(probitsim$a4); probitsim$a4 = NULL

# outcome variable
wgt3 = probitsim$wgt3; probitsim$wgt3 = NULL

# data transformation
probitsim$cage<-probitsim$age-mean(probitsim$age)
probitsim$cage2<-probitsim$cage^2
probitsim$cwgt0<-probitsim$wgt0-mean(probitsim$wgt0)
probitsim$cwgt02<-probitsim$cwgt0^2

# unused variable
probitsim$id = NULL


#### estimators ####

normalize = function(W, A) W / ifelse(A, sum(W[A]), sum(W[!A]))

get_output = function(W.ATE, W.ATT, X, Y, y0, y1, Z,
                           mask.ATE = rep(T, length(Y)),
                           mask.ATT = rep(T, length(Y))) {
  if (!all(sapply(list(Z, mask.ATE, mask.ATT), is.logical)))
    stop(paste(paste(
      c("'Z'", "'mask.ATE'", "'mask.ATT'")[sapply(list(Z, mask.ATE, mask.ATT), is.logical)],
      collapse = ", "), "must be logical"))
  
  ATE.linreg = try(lm(Y ~ Z, weights = W.ATE*mask.ATE)$coefficients[[2]])
  ATT.linreg = try(lm(Y ~ Z, weights = W.ATT*mask.ATT)$coefficients[[2]])
  
  return(
    c("ATE" = sum((2*Z-1)*W.ATE*Y*mask.ATE),
      "ATE.dr" = mean((y1-y0)[mask.ATE]) + sum(mask.ATE*W.ATE*ifelse(Z,(Y-y1),-(Y-y0))),
      "ATE.linreg" = if (inherits(ATE.linreg, "try-error")) NA else ATE.linreg,
      "ATT" = mean(Y[Z & mask.ATT]) - sum(W.ATT*Y*mask.ATT*!Z),
      "ATT.dr" = mean((Y-y0)[Z & mask.ATT]) - sum(W.ATT*(Y-y0)*mask.ATT*!Z),
      "ATE.linreg" = if (inherits(ATT.linreg, "try-error")) NA else ATT.linreg
    )
  )
}

#### comparison ####


# import our weights
load("output.Rdata")


plot_density = function(W, A, main = "", prob = F, add = F, legend = T) {
  d1 = density(W[A]); if(prob) d1$y = d1$y / sum(d1$y) * 512 / (diff(range(d1$x)))
  d0 = density(W[!A]); if(prob) d0$y = d0$y / sum(d0$y) * 512 / (diff(range(d0$x)))
  (if(add) lines else plot)(d1, xlim=range(c(d0$x,d1$x)), ylim=range(c(d0$y,d1$y)), col="blue",
                            main = main, lwd = 2)
  lines(d0, col="red", lwd=2)
  if (legend) legend("topright", c("control","treated"), col=c("red", "blue"), lty = 1, lwd=2)
}



##### response surface #####

###### log reg ######

logreg_rs_train = function(X, Y, Z) {
  list(
    "y0" = lm(
      Y ~ factor(educ) + factor(smoke) + factor(allergy) + cage + cage2 + factor(location),
      data = X,
      subset = !Z
    ),
    "y1" = lm(
      Y ~ factor(educ) + factor(smoke) + factor(allergy) + cage + cage2 + factor(location),
      data = X,
      subset = Z
    )
  )
}

logreg_rs_pred = function(X, trained.model) {
  lapply(trained.model,
         function(mdl) predict(mdl, X, type = "response") |> setNames(NULL)
  )
}

###### rf ######

rf_rs_train = function(X, Y, Z) {
  trctrl = trainControl(method="cv", number=5)
  
  list(
    "y0" = train(
      Y ~ factor(educ) + factor(smoke) + factor(allergy) + cage + cage2 + factor(location),
      data = data.frame(X, Y = Y), subset = !Z, method = "rf", metric = "RMSE", trControl = trctrl
    ),
    "y1" = train(
      Y ~ factor(educ) + factor(smoke) + factor(allergy) + cage + cage2 + factor(location),
      data = data.frame(X, Y = Y), subset = Z, method = "rf", metric = "RMSE", trControl = trctrl
    )
  )
}

rf_rs_pred = function(X, trained.model) {
  lapply(trained.model, function(m) predict(m, newdata = X) |> setNames(NULL))
}

##### calcul y0, y1 #####
# rs = list()
# rs[[scenario]]$logreg = logreg_rs_train(X, Y, A) |> logreg_rs_pred(X=X)
# rs[[scenario]]$rf = rf_rs_train(X, Y, A) |>     rf_rs_pred(X=X)
# save(rs, file = "rs.Rdata")
load("rs.Rdata")

#### Compute Stat ####

par(mfrow = c(2,3))
output = lapply(
  c("a2", "a3 | a1=0", "a3 | a1=1"),
  function(scenario) {
    mask = switch(scenario,
                  "a2" = rep(T, nrow(probitsim)),
                  "a3 | a1=0" = !a1,
                  "a3 | a1=1" = a1)
    X = probitsim[mask, ]
    Y = wgt3[mask]
    A = switch(scenario, "a2" = a2, a3)[mask]
    
    y0_logreg = rs[[scenario]]$logreg$y0
    y1_logreg = rs[[scenario]]$logreg$y1
    y0_rf = rs[[scenario]]$rf$y0
    y1_rf = rs[[scenario]]$rf$y1
    
    out = list()
    ##### IPTW #####
    ###### log reg ######
    W = data[[scenario]]$`AIPW_log reg`
    plot_density(data[[scenario]]$`AIPW_log reg`$ATE, A,
                 main = "IPTW log reg", legend=F)
    
    # no post processing IPTW
    out[["IPTW_logreg"]] = get_output(
      W.ATE = W$ATE,
      W.ATT = W$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_logreg, y1 = y1_logreg
    )
    names(out[["IPTW_logreg"]])[c(2,5)] = c("ATE.dr.logreg", "ATT.dr.logreg")
    out[["IPTW_logreg"]] = append(out[["IPTW_logreg"]], c("ATE.dr.rf"=NA), after=2)
    out[["IPTW_logreg"]] = append(out[["IPTW_logreg"]], c("ATT.dr.rf"=NA), after=6)
    
    # IPTW w/ norm
    W.ATE.norm = normalize(W$ATE, A)
    W.ATT.norm = normalize(W$ATT, A)
    out[["IPTW_logreg.norm"]] = get_output(
      W.ATE = W.ATE.norm,
      W.ATT = W.ATT.norm,
      X = X, Y = Y, Z = A,
      y0 = y0_logreg, y1 = y1_logreg
    )
    names(out[["IPTW_logreg.norm"]])[c(2,5)] = c("ATE.dr.logreg", "ATT.dr.logreg")
    out[["IPTW_logreg.norm"]] = append(out[["IPTW_logreg.norm"]], c("ATE.dr.rf"=NA), after=2)
    out[["IPTW_logreg.norm"]] = append(out[["IPTW_logreg.norm"]], c("ATT.dr.rf"=NA), after=6)
    
    # IPTW w/ trim
    idx.trim.ATE = W$ATE < quantile(W$ATE,0.99)
    idx.trim.ATT = (W$ATT*!A) < quantile(W$ATT[!A], .99)
    
    out[["IPTW_logreg.trim"]] = get_output(
      W.ATE = W$ATE,
      W.ATT = W$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_logreg, y1 = y1_logreg,
      mask.ATE = idx.trim.ATE,
      mask.ATT = idx.trim.ATT
    )
    names(out[["IPTW_logreg.trim"]])[c(2,5)] = c("ATE.dr.logreg", "ATT.dr.logreg")
    out[["IPTW_logreg.trim"]] = append(out[["IPTW_logreg.trim"]], c("ATE.dr.rf"=NA), after=2)
    out[["IPTW_logreg.trim"]] = append(out[["IPTW_logreg.trim"]], c("ATT.dr.rf"=NA), after=6)
    
    #IPTW w/ norm then trim
    out[["IPTW_logreg.norm.trim"]] = get_output(
      W.ATE = W.ATE.norm,
      W.ATT = W.ATT.norm,
      X = X, Y = Y, Z = A,
      y0 = y0_logreg, y1 = y1_logreg,
      mask.ATE = W.ATE.norm < quantile(W.ATE.norm, .99),
      mask.ATT = (W.ATT.norm*!A) < quantile(W.ATT.norm[!A], .99)
    )
    names(out[["IPTW_logreg.norm.trim"]])[c(2,5)] = c("ATE.dr.logreg", "ATT.dr.logreg")
    out[["IPTW_logreg.norm.trim"]] = append(out[["IPTW_logreg.norm.trim"]], c("ATE.dr.rf"=NA), after=2)
    out[["IPTW_logreg.norm.trim"]] = append(out[["IPTW_logreg.norm.trim"]], c("ATT.dr.rf"=NA), after=6)
    
    #IPTW w/ trim then norm
    out[["IPTW_logreg.trim.norm"]] = get_output(
      W.ATE = normalize(W$ATE*idx.trim.ATE, A),
      W.ATT = normalize(W$ATT * idx.trim.ATT, A),
      X = X, Y = Y, Z = A,
      y0 = y0_logreg, y1 = y1_logreg,
      mask.ATE = idx.trim.ATE,
      mask.ATT = idx.trim.ATT
    )
    names(out[["IPTW_logreg.trim.norm"]])[c(2,5)] = c("ATE.dr.logreg", "ATT.dr.logreg")
    out[["IPTW_logreg.trim.norm"]] = append(out[["IPTW_logreg.trim.norm"]], c("ATE.dr.rf"=NA), after=2)
    out[["IPTW_logreg.trim.norm"]] = append(out[["IPTW_logreg.trim.norm"]], c("ATT.dr.rf"=NA), after=6)
    
    
    ###### rf ######
    
    W = data[[scenario]]$`AIPW_rf`
    plot_density(data[[scenario]]$`AIPW_rf`$ATE, A,
                 main = "IPTW rf", legend=F)
    
    # no post processing IPTW
    out[["IPTW_rf"]] = get_output(
      W.ATE = W$ATE,
      W.ATT = W$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_rf, y1 = y1_rf
    )
    names(out[["IPTW_rf"]])[c(2,5)] = c("ATE.dr.rf", "ATT.dr.rf")
    out[["IPTW_rf"]] = append(out[["IPTW_rf"]], c("ATE.dr.logreg"=NA), after=1)
    out[["IPTW_rf"]] = append(out[["IPTW_rf"]], c("ATT.dr.logreg"=NA), after=5)
    
    # IPTW w/ norm
    W.ATE.norm = normalize(W$ATE, A)
    W.ATT.norm = normalize(W$ATT, A)
    out[["IPTW_rf.norm"]] = get_output(
      W.ATE = W.ATE.norm,
      W.ATT = W.ATT.norm,
      X = X, Y = Y, Z = A,
      y0 = y0_rf, y1 = y1_rf
    )
    names(out[["IPTW_rf.norm"]])[c(2,5)] = c("ATE.dr.rf", "ATT.dr.rf")
    out[["IPTW_rf.norm"]] = append(out[["IPTW_rf.norm"]], c("ATE.dr.logreg"=NA), after=1)
    out[["IPTW_rf.norm"]] = append(out[["IPTW_rf.norm"]], c("ATT.dr.logreg"=NA), after=5)
    
    # IPTW w/ trim
    idx.trim.ATE = W$ATE < quantile(W$ATE,0.99)
    idx.trim.ATT = (W$ATT*!A) < quantile(W$ATT[!A], .99)
    
    out[["IPTW_rf.trim"]] = get_output(
      W.ATE = W$ATE,
      W.ATT = W$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_rf, y1 = y1_rf,
      mask.ATE = idx.trim.ATE,
      mask.ATT = idx.trim.ATT
    )
    names(out[["IPTW_rf.trim"]])[c(2,5)] = c("ATE.dr.rf", "ATT.dr.rf")
    out[["IPTW_rf.trim"]] = append(out[["IPTW_rf.trim"]], c("ATE.dr.logreg"=NA), after=1)
    out[["IPTW_rf.trim"]] = append(out[["IPTW_rf.trim"]], c("ATT.dr.logreg"=NA), after=5)
    
    #IPTW w/ norm then trim
    out[["IPTW_rf.norm.trim"]] = get_output(
      W.ATE = W.ATE.norm,
      W.ATT = W.ATT.norm,
      X = X, Y = Y, Z = A,
      y0 = y0_rf, y1 = y1_rf,
      mask.ATE = W.ATE.norm < quantile(W.ATE.norm, .99),
      mask.ATT = (W.ATT.norm*!A) < quantile(W.ATT.norm[!A], .99)
    )
    names(out[["IPTW_rf.norm.trim"]])[c(2,5)] = c("ATE.dr.rf", "ATT.dr.rf")
    out[["IPTW_rf.norm.trim"]] = append(out[["IPTW_rf.norm.trim"]], c("ATE.dr.logreg"=NA), after=1)
    out[["IPTW_rf.norm.trim"]] = append(out[["IPTW_rf.norm.trim"]], c("ATT.dr.logreg"=NA), after=5)
    
    #IPTW w/ trim then norm
    out[["IPTW_rf.trim.norm"]] = get_output(
      W.ATE = normalize(W$ATE*idx.trim.ATE, A),
      W.ATT = normalize(W$ATT * idx.trim.ATT, A),
      X = X, Y = Y, Z = A,
      y0 = y0_rf, y1 = y1_rf,
      mask.ATE = idx.trim.ATE,
      mask.ATT = idx.trim.ATT
    )
    names(out[["IPTW_rf.trim.norm"]])[c(2,5)] = c("ATE.dr.rf", "ATT.dr.rf")
    out[["IPTW_rf.trim.norm"]] = append(out[["IPTW_rf.trim.norm"]], c("ATE.dr.logreg"=NA), after=1)
    out[["IPTW_rf.trim.norm"]] = append(out[["IPTW_rf.trim.norm"]], c("ATT.dr.logreg"=NA), after=5)
    
    
    ##### TLF weights #####
    plot_density(data[[scenario]]$`TLF_log reg`$ATE, A,
                 main = "TLF", legend=F)
    
    out[["TLF"]] = get_output(
      W.ATE = data[[scenario]]$`TLF_log reg`$ATE,
      W.ATT = data[[scenario]]$`TLF_log reg`$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_logreg, y1 = y1_logreg
    )
    names(out[["TLF"]])[c(2,5)] = c("ATE.dr.logreg", "ATT.dr.logreg")
    out[["TLF"]] = append(out[["TLF"]], c("ATE.dr.rf"=NA), after=2)
    out[["TLF"]] = append(out[["TLF"]], c("ATT.dr.rf"=NA), after=6)
    
    out[["TLF"]][c(3,7)] = get_output(
      W.ATE = data[[scenario]]$`TLF_rf`$ATE,
      W.ATT = data[[scenario]]$`TLF_rf`$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_rf, y1 = y1_rf
    )[c(2,5)]
    
    ##### EB weights #####
    plot_density(data[[scenario]]$`EB_log reg`$ATE, A,
                 main = "EB", legend=F)
    
    out[["EB"]] = get_output(
      W.ATE = data[[scenario]]$`EB_log reg`$ATE,
      W.ATT = data[[scenario]]$`EB_log reg`$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_logreg, y1 = y1_logreg
    )
    names(out[["EB"]])[c(2,5)] = c("ATE.dr.logreg", "ATT.dr.logreg")
    out[["EB"]] = append(out[["EB"]], c("ATE.dr.rf"=NA), after=2)
    out[["EB"]] = append(out[["EB"]], c("ATT.dr.rf"=NA), after=6)
    
    out[["EB"]][c(3,7)] = get_output(
      W.ATE = data[[scenario]]$`EB_rf`$ATE,
      W.ATT = data[[scenario]]$`EB_rf`$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_rf, y1 = y1_rf
    )[c(2,5)]
    
    ##### KOM weights #####
    plot_density(data[[scenario]]$`KOM_log reg`$ATE, A,
                 main = "KOM", legend=F)
    
    out[["KOM"]] = get_output(
      W.ATE = data[[scenario]]$`KOM_log reg`$ATE,
      W.ATT = data[[scenario]]$`KOM_log reg`$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_logreg, y1 = y1_logreg
    )
    names(out[["KOM"]])[c(2,5)] = c("ATE.dr.logreg", "ATT.dr.logreg")
    out[["KOM"]] = append(out[["KOM"]], c("ATE.dr.rf"=NA), after=2)
    out[["KOM"]] = append(out[["KOM"]], c("ATT.dr.rf"=NA), after=6)
    
    out[["KOM"]][c(3,7)] = get_output(
      W.ATE = data[[scenario]]$`KOM_rf`$ATE,
      W.ATT = data[[scenario]]$`KOM_rf`$ATT,
      X = X, Y = Y, Z = A,
      y0 = y0_rf, y1 = y1_rf
    )[c(2,5)]
    
    # create custom legend
    plot.new()
    lines(c(0,.2),c(.66,.66), col='red', lwd=2)
    text(.3,.66,"control group", adj=0)
    lines(c(0,.2),c(.33,.33), col='blue', lwd=2)
    text(.3,.33,"treated group", adj=0)
    mtext(paste("Scenario :", scenario), side = 3, line = - 1.5, outer = TRUE)
    return(do.call(rbind, out) |> round(0))
  }
)
output[[1]]

