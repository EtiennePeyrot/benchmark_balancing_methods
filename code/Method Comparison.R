# Study estimation ATE / ATT

#### Initialization ####

# setwd("G:/Mon Drive/Study/PhD/Article 1 - Benchmarking Balancing Method/Code/only server")
rm(list = ls())
graphics.off()
set.seed(123)

# File import

#install.packages(c("magrittr", "randomForest", "caret", "WeightIt", "osqp", "sandwich",
#                   "Mestim", "mvtnorm", "doParallel"))

source("Method Definition.R")
source("Crossfit and Compute Method.R")

# library import

library(mvtnorm) # used in fct : GenerateData
library(doParallel) # multi processing

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

GenerateData = function(N, sample_size, trt_rarity, confdg_lvl) {
  # Generate N samples of size sample_size with various level of confounding and
  # treatment rarity.
  #
  # args :
  #       - N : numeric, the number of sample to generate
  #       - sample_size : integer, the size of each sample
  #       - trt_rarity : str, the rarity of the treatement, can be
  #                      "common" (~35%), "rare" (~15%), "very rare" (~5%),
  #       - confdg : str, the level of confounding for the simulation. can be
  #                  "low", "moderate" or "high"
  
  if (N <= 0 && N != floor(N))
    stop("N must be a positive integer")
  
  if (sample_size <= 0 && sample_size != floor(sample_size))
    stop("sample_size must be a positive integer")
  
  if (!(trt_rarity %in% c("common", "rare", "very rare")))
    stop("trt_rarity must be 'common', 'rare' or 'very rare'")
  
  if (!(confdg_lvl %in% c("low", "moderate", "high")))
    stop("confdg must be 'low', 'moderate' or 'high'")
  
  # covariates
  
  Sigma = diag(rep(1,10)) # covariance matrix
  Sigma[5,1] = Sigma[1,5] = 0.2
  Sigma[8,3] = Sigma[3,8] = 0.2
  Sigma[6,2] = Sigma[2,6] = 0.9
  Sigma[9,4] = Sigma[4,9] = 0.9
  
  X = rmvnorm(N*sample_size, rep(0, 10), Sigma)
  X[ ,c(1,3,5,6,8,9)] = as.numeric(X[ ,c(1,3,5,6,8,9)] > 0)
  
  
  # treatment assignment
  
  b0 = switch(trt_rarity, "common" = -1.84, "rare" = -4.12, "very rare" = -6.5)
  b = c(0.8, -0.25, 0.6, -0.4, -0.8, -0.5, 0.7, 0, 0, 0)
  A = rbinom(N*sample_size, 1, expit( b0 + 2.25 * ( X %*% b + 0.5*X[,1]*X[,2]^2) ))
  
  
  # Outcome
  
  g = switch(confdg_lvl, "low" = 1, "moderate" = 2.25, "high" = 5)
  a0 = switch(confdg_lvl, "low" = -1.5, "moderate" = -2.22, "high" = -4.1)
  a = c(0.3, -0.36, -0.73, -0.2, 0, 0, 0, 0.71, -0.19, 0.26)
  Y = rbinom(N * sample_size, 1, expit(a0 + g * ( X %*% a + 0.5*X[,3]*X[,4]^2) ))
  
  X = data.frame(X)
  X[,-c(1,3,5,6,8,9)] = scale(X[,-c(1,3,5,6,8,9)])
  
  output = list("X" = X,
                "Z" = A,
                "Y" = Y,
                "confdg_lvl" = confdg_lvl,
                "trt_rarity" = trt_rarity)
  
  return(output)
}

#### Nuisance Functions ####

# propensity score 

## random forest

ps_rf_train = function(X, Z) Propensity_Score(
  X = X, Z = Z, method = "rf", target = "train", trained.model = NULL)

ps_rf_pred = function(X, Z, trained.model) Propensity_Score(
  X = X, Z = Z, method = "rf", target = "predict", trained.model = trained.model)

## well specified logistic regression

ps_true_train = function(X, Z) Propensity_Score(
  X = X, Z = Z, method = "true", target = "train", trained.model = NULL)

ps_true_pred = function(X, Z, trained.model) Propensity_Score(
  X = X, Z = Z, method = "true", target = "predict",trained.model = trained.model)

## misspecified logistic regression

ps_misspecified_train = function(X, Z) Propensity_Score(
  X = X, Z = Z, method = "misspecified", target = "train", trained.model = NULL)

ps_misspecified_pred = function(X, Z, trained.model) Propensity_Score(
  X = X, Z = Z, method = "misspecified", target = "predict", trained.model = trained.model)


#response surfaces

## random forest

rs_rf_train = function(X, Y, Z) Response_Surface(
  X = X, Y = Y, Z = Z, method = "rf", target = "train", trained.model = NULL)

rs_rf_pred = function(X, Y, Z, trained.model) Response_Surface(
  X = X, Y = Y, Z = Z, method = "rf", target = "predict", trained.model = trained.model)

## well specified logistic regression

rs_true_train = function(X, Y, Z) Response_Surface(
  X = X, Y = Y, Z = Z, method = "true", target = "train", trained.model = NULL)

rs_true_pred = function(X, Y, Z, trained.model) Response_Surface(
  X = X, Y = Y, Z = Z, method = "true", target = "predict", trained.model = trained.model)

## misspecified logistic regression

rs_misspecified_train = function(X, Y, Z) Response_Surface(
  X = X, Y = Y, Z = Z, method = "misspecified", target = "train", trained.model = NULL)

rs_misspecified_pred = function(X, Y, Z, trained.model) Response_Surface(
  X = X, Y = Y, Z = Z, method = "misspecified", target = "predict", trained.model = trained.model)


####     Parameters     ####

# Data Generation
N = 5000 # number of repetition : 7300 ou +
SAMPLE_SIZE = c(2000) # c(250,500,1000,2000)
TRT_RARITY = c("very rare") # c("common", "rare", "very rare")
CONFDG_LVL = c("moderate", "high") # c("low", "moderate", "high")

# Methods
METHODS = c("AIPW", "EB", "KOM", "TLF") # c("AIPW", "EB", "KOM", "TLF")
MODELS = c("rf") # c("true", "misspecified", "rf")
crossfit_methods = FALSE # sample_size >= 1000
nbr_crossfit = 30
nbr_split = 5

# Parallelization
parallel = TRUE
ncores = 20 # detectCores(logical = TRUE)
external_script = c("Method Definition.R", "Crossfit and Compute Method.R")

# Other
save_output = TRUE
generate_Data = FALSE
dataset_per_file = 50


#### Data Generation ####

if (generate_Data) {
  set.seed(123)
  cat("Data generation...\n")
  nbr_file = length(SAMPLE_SIZE) * length(TRT_RARITY) * length(CONFDG_LVL) * ceiling(N/dataset_per_file)
  cptr_file = 0
  for (sample_size in SAMPLE_SIZE) {
    for (trt_rarity in TRT_RARITY) {
      for (confdg_lvl in CONFDG_LVL) {
        for (i in 1:ceiling(N/dataset_per_file)) {
          cptr_file = cptr_file + 1
          cat("\r",format(cptr_file,width = ceiling(log10(nbr_file)),justify = "right"),
              " file out of ",nbr_file," created",sep = "")
          filename = paste(sample_size, "obs", trt_rarity,
                           "trt", confdg_lvl, "confdg",i,
                           sep = "_")
          
          sim = GenerateData(N = dataset_per_file,
                             sample_size = sample_size,
                             trt_rarity = trt_rarity,
                             confdg_lvl = confdg_lvl)
          
          save(sim, file = paste0("Data/",filename))
        }
      }
    }
  }
  cat("\nDone!")
}


#### Balancing Methods ####

Balancing_Method = list()
for (method in METHODS) {
  for (model in MODELS) {
    for (cf in if (crossfit_methods) c(T, F) else c(F)) {
      Balancing_Method[[paste0(method, ifelse(cf, "_cf_", "_"), model)]] =
        list(
          # eval(substitute(function(i) f(i, bool), list(bool = b)))
          "method" = switch(method,
                             "AIPW" = eval(substitute(
                               function(X, Y, Z, propensity_score,response_surface)
                                 AIPW(X, Y, Z, propensity_score, response_surface, var = bool),
                               list(bool = !cf))),
                             
                             "EB"   = eval(substitute(
                               function(X,Y,Z,response_surface)
                                 EnergyB(X,Y,Z,response_surface, var = bool),
                               list(bool = !cf))),
                             
                             "KOM"  = eval(substitute(
                               function(X,Y,Z,response_surface)
                                 KOM(X,Y,Z,response_surface, var = bool),
                               list(bool = !cf))),
                             
                             "TLF"  = eval(substitute(
                               function(X,Y,Z,response_surface) Tailored_loss_function(
                                 X,Y,Z,response_surface,100,.001,NULL, var = bool),
                               list(bool = !cf)))),
          
          "nuisance_fct" = if (method != "AIPW") list(
            "response_surface" = list(
              "train" = switch(model,
                               "true"         = rs_true_train,
                               "misspecified" = rs_misspecified_train,
                               "rf"           = rs_rf_train),
              "predict" = switch(model,
                                 "true"         = rs_true_pred,
                                 "misspecified" = rs_misspecified_pred,
                                 "rf"           = rs_rf_pred)
            )
          ) else  list(
            "response_surface" = list(
              "train" = switch(model,
                               "true"         = rs_true_train,
                               "misspecified" = rs_misspecified_train,
                               "rf"           = rs_rf_train),
              "predict" = switch(model,
                                 "true"         = rs_true_pred,
                                 "misspecified" = rs_misspecified_pred,
                                 "rf"           = rs_rf_pred)
            ),
            "propensity_score" = list(
              "train" = switch(model,
                               "true"         = ps_true_train,
                               "misspecified" = ps_misspecified_train,
                               "rf"           = ps_rf_train),
              "predict" = switch(model,
                                 "true"         = ps_true_pred,
                                 "misspecified" = ps_misspecified_pred,
                                 "rf"           = ps_rf_pred)
            )
          ),
          
          "crossfit" = cf,
          "nbr_crossfit" = nbr_crossfit,
          "nbr_split" = nbr_split
        )
    }
  }
}

error_msg = check_method(Balancing_Method)
if (error_msg != "")
  stop(error_msg)


#### Core Code ####

for (sample_size in SAMPLE_SIZE) {
  for (trt_rarity in TRT_RARITY) {
    for (confdg_lvl in CONFDG_LVL) {
      for (split in if (confdg_lvl == "moderate") 21:100 else 1:100) { # 1:((N - 1) %/% dataset_per_file + 1)) { # 1:(N %/% dataset_per_file + 1)
        
	# if (!inherits(try(load(paste0("Output/", paste(sample_size, "obs", confdg_lvl, "confounding", trt_rarity, "trt", split)))),"try-error")) next
        # load data
        filename = paste(sample_size, "obs", trt_rarity, "trt",
                         confdg_lvl, "confdg", split, sep = "_")
        load(paste0("Data/", filename))
        
        if (!parallel) {
          df = lapply(
            1:dataset_per_file,
            function(i) {
              sample_index = ((i - 1) %% dataset_per_file) * sample_size + 1:sample_size
              
              # compute the methods
              compute_method(
                X = sim$X[sample_index, ],
                Y = sim$Y[sample_index],
                Z = sim$Z[sample_index],
                Balancing_Method = Balancing_Method)
            }) %>% {
              # concatenate and rename results
              name = names(.[[1]])
              out = lapply(name, function(n) lapply(., `[[`, n) %>%
                             do.call(rbind, .))
              `names<-`(out, name)
            }
        } else {
          # create cluster
          cl = makeCluster(ncores)
          registerDoParallel(cl)
          
          df = foreach(i = 1:dataset_per_file,
                       .final = .mycombine,
                       .packages = c("magrittr", "caret", "WeightIt", "osqp", "sandwich", "Mestim")
          ) %dopar% {
            
            # import script
            my_source(external_script)
            
            # compute the methods
            sample_index = ((i - 1) %% dataset_per_file) * sample_size + 1:sample_size
            compute_method(
              X = sim$X[sample_index, ],
              Y = sim$Y[sample_index],
              Z = sim$Z[sample_index],
              Balancing_Method = Balancing_Method)
          }
          # close clusters
          stopCluster(cl)
        }
        
        # save results in 'Data'
        Data = list('X' = sim$X,'Y' = sim$Y, 'Z' = sim$Z, 'df' = df,
                    'confdg_lvl' = confdg_lvl, 'trt_rarity' = trt_rarity)
        
        if (save_output) paste0("Output/", sample_size, " obs ",
                              confdg_lvl, " confounding ",
                              trt_rarity, " trt ", split) %>%
          save(Data, file = .)
      }
    }
  }
}