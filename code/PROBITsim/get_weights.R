#### initialization ####

setwd("C:/Users/skoua/project benchmark balancing/PROBITsim")
rm(list = ls())
graphics.off()
set.seed(123)

# file import

source("Method Definition.R")
source("Crossfit and Compute Method.R")

# library import

library(foreign) # import dta file
library(randomForest) # random forest

# data import

probitsim = as.data.frame(read.dta("PROBITsim2018_v12.dta", convert.factors = FALSE))

#### Functions ####

my_source = function(scripts) for (script in c(scripts)) source(script)


expit = function(X) 1 / (1 + exp(-X)) # sigmoid function

.mycombine = function(l) {
  name = names(l[[1]])
  out = lapply(name,
               function(n) do.call(rbind, lapply(l, `[[`, n)) )
  names(out) = name
  return(out)
}


#### pre processing ####

# treatment variables
probitsim$a1 = as.logical(probitsim$a1)
probitsim$a2 = as.logical(probitsim$a2)
probitsim$a3 = as.logical(probitsim$a3)
probitsim$a4 = as.logical(probitsim$a4)

# outcome variable
wgt3 = probitsim$wgt3; probitsim$wgt3 = NULL

# data transformation
probitsim$cage<-probitsim$age-mean(probitsim$age)
probitsim$cage2<-probitsim$cage^2
probitsim$cwgt0<-probitsim$wgt0-mean(probitsim$wgt0)
probitsim$cwgt02<-probitsim$cwgt0^2

# unused variable
probitsim$id = NULL


#### nuissance functions ####

##### propensity score #####

###### log reg ######

logreg_ps_a2_train = function(X, Z) {
  glm(
    Z ~ factor(location) + factor(educ) + cage + cage2 + factor(smoke) + factor(allergy),
    family = binomial,
    data = X
  )
}

logreg_ps_a3_train = function(X, Z) {
  glm(
    Z ~ a2 + factor(location) + factor(educ) + factor(smoke) + factor(allergy) +
      factor(cesarean) + factor(sex) + cage + cage2 + cwgt0 + cwgt02,
    family = binomial,
    data = X
  )
}

logreg_ps_pred = function(X, trained.model) {
  predict(trained.model, newdata = X, type = "response") |> setNames(NULL)
}

###### rf ######

rf_ps_a2_train = function(X, Z) {
  trctrl = trainControl(method="cv", number=5)
  train(
    Z ~ factor(location) + factor(educ) + cage + cage2 + factor(smoke) + factor(allergy),
    data = data.frame(X, Z = as.numeric(Z)), method = "rf", metric = "RMSE", trControl = trctrl
  )
}

rf_ps_a3_train = function(X, Z) {
  trctrl = trainControl(method="cv", number=5)
  train(
    Z ~ a2 + factor(location) + factor(educ) + factor(smoke) + factor(allergy) +
      factor(cesarean) + factor(sex) + cage + cage2 + cwgt0 + cwgt02,
    data = data.frame(X, Z = as.numeric(Z)), method = "rf", metric = "RMSE", trControl = trctrl
  )
}

rf_ps_pred = function(X, trained.model) 
  predict(trained.model, newdata = X, type = "raw") |> setNames(NULL)



##### response surface #####

###### log reg ######

logreg_rs_a2_train = function(X, Y, Z) {
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

logreg_rs_a3_train = function(X, Y, Z) {
  list(
    "y0" = lm(
      
      Y ~ a2 + factor(location) + factor(educ) + factor(smoke) + factor(allergy) +
        factor(cesarean) + factor(sex) + cage + cage2 + cwgt0 + cwgt02,
      data = X,
      subset = !Z
    ),
    "y1" = lm(
      Y ~ a2 + factor(location) + factor(educ) + factor(smoke) + factor(allergy) +
        factor(cesarean) + factor(sex) + cage + cage2 + cwgt0 + cwgt02,
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

rf_rs_a2_train = function(X, Y, Z) {
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

rf_rs_a3_train = function(X, Y, Z) {
  trctrl = trainControl(method="cv", number=5)
  
  list(
    "y0" = train(
      Y ~ a2 + factor(location) + factor(educ) + factor(smoke) + factor(allergy) +
        factor(cesarean) + factor(sex) + cage + cage2 + cwgt0 + cwgt02,
      data = data.frame(X, Y = Y), subset = !Z, method = "rf", metric = "RMSE", trControl = trctrl
    ),
    "y1" = train(
      Y ~ a2 + factor(location) + factor(educ) + factor(smoke) + factor(allergy) +
        factor(cesarean) + factor(sex) + cage + cage2 + cwgt0 + cwgt02,
      data = data.frame(X, Y = Y), subset = Z, method = "rf", metric = "RMSE", trControl = trctrl
    )
  )
}

rf_rs_pred = function(X, trained.model) {
  lapply(trained.model, function(m) predict(m, newdata = X) |> setNames(NULL))
}

#### create estimation methods ####

Balancing_Method = list()
for (scenario in c("a2", "a3")) {
  for (bal_met in c("AIPW", "EB", "KOM", "TLF")) {
    for (setting in c("log reg", "rf")) {
      Balancing_Method[[scenario]][[paste0(bal_met, "_", setting)]] =
        list(
          "method" = switch(bal_met,
                            "AIPW" =
                              function(X, Y, Z, propensity_score,response_surface)
                                AIPW(X, Y, Z, propensity_score, response_surface, var = T),
                            
                            "EB"   =
                              function(X,Y,Z,response_surface)
                                EnergyB(X,Y,Z,response_surface, var = T),
                            
                            "KOM"  =
                              function(X,Y,Z,response_surface)
                                KOM(X,Y,Z,response_surface, var = T),
                            
                            "TLF"  =
                              function(X,Y,Z,response_surface) Tailored_loss_function(
                                X,Y,Z,response_surface,nbr_fold=5, var = T)
                              ),
          
          "nuisance_fct" = c(list(
            "response_surface" = list(
              "train" = switch(setting,
                               "log reg" = switch(scenario, "a2" = logreg_rs_a2_train, "a3" = logreg_rs_a3_train),
                               "rf"      = switch(scenario, "a2" = rf_rs_a2_train, "a3" = rf_rs_a3_train)),
              "predict" = switch(setting,
                                 "log reg" = logreg_rs_pred,
                                 "rf"      = rf_rs_pred)
            )),
            if (bal_met != "AIPW") NULL else list(
              "propensity_score" = list(
                "train" = switch(setting,
                                 "log reg" = switch(scenario, "a2" = logreg_ps_a2_train, "a3" = logreg_ps_a3_train),
                                 "rf"      = switch(scenario, "a2" = rf_ps_a2_train, "a3" = rf_ps_a3_train)),
                "predict" = switch(setting,
                                   "log reg" = logreg_ps_pred,
                                   "rf"      = rf_ps_pred)
                )
              )
            ),
          "crossfit" = F
        )
    }
  }
}

#### estimate the treatment effect of a2
res1=compute_method(probitsim, wgt3, probitsim$a2, Balancing_Method[["a2"]])
for (name in names(res1)) res1[[name]]$name = name
cat("res1 done!")

#### estimate the treatment effect of a3 under a1=0
filter = !probitsim$a1
res2=compute_method(probitsim[filter, ], wgt3[filter], probitsim$a3[filter], Balancing_Method[["a3"]])
for (name in names(res2)) res2[[name]]$name = name
cat("res2 done!")

#### estimate the treatment effect of a3 under a1=1
filter = probitsim$a1
res2=compute_method(probitsim[filter, ], wgt3[filter], probitsim$a3[filter], Balancing_Method[["a3"]])
for (name in names(res3)) res3[[name]]$name = name
cat("res3 done!")

# save outputs
data = list("a2" = res1,
            "a3 | a1=0" = res2,
            "a3 | a1=1" = res3)

save(data, file = "output.Rdata")
cat("job done!")
