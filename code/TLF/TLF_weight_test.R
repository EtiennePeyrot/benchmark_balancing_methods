rm(list = ls())
graphics.off()
# set.seed(123)


# import
library(ggplot2)
source("confidence_interval.R")
source("TLF_weight.R")
source("../GenerateData.R")


# function
norm = function(w,a) w / ifelse(a, sum(w[a]), sum(w[!a]))

# input
n = 2000 # sample size
prop_trt = "moderate" # c("low", "moderate", "high"),
confdg_lvl = "high" # c("low", "moderate", "high")

# generate data
data = GenerateData(n = n, prop_trt = prop_trt, confdg_lvl = confdg_lvl)
X = data$X
A = as.logical(data$A)
A_ = ifelse(A,"treated","control")
Y = data$Y

# pre computed hyper parameters
lambda = array(data = c(0.01 , 0.001, 0.01 , 0.001, 0.001, 1e-04, 0.001, 1e-04, 0.01 ,
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
lambda = as.list(lambda[,as.character(n), prop_trt, confdg_lvl])


# compute weights
time = system.time({
  D1 = as.matrix(dist(X, method = "manhattan", diag = T, upper = T))
  sigma = 1 / median(D1)
  K.eigen = eigen(exp(-sigma * D1))
  res = TLF_weight(X = X, A = A, sigma = sigma, lambda = lambda,
                 K.eigen = K.eigen, intercept = T)
  
  res_ate = res$model$ATE
  res_att = res$model$ATT
  # Propensity block for variance sandwich estimator
  ps_block_ate = propensity_blocks_cbpstlf(A = A, U = res_ate$U, eta = res_ate$eta,
                                           d = res_ate$d, lambda = lambda$ATE,
                                           target="ATE")
  ps_block_att = propensity_blocks_cbpstlf(A = A, U = res_att$U, eta = res_att$eta,
                                           d = res_att$d, lambda = lambda$ATT,
                                           target="ATT")
  
  wls_ate = wls_ci(Y, A, ps_block = ps_block_ate, target = "ATE", alpha = 0.05)
  wls_att = wls_ci(Y, A, ps_block = ps_block_att, target = "ATT", alpha = 0.05)
  
  # Prognostic models
  outcome_block = outcome_blocks_binomial(Y = Y, A = A, X = X)
  
  # AIPW
  dr_ate = dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_ate, target = "ATE", alpha = 0.05)
  dr_att = dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block_att, target = "ATT", alpha = 0.05)
  
  estimate = list(
    WLS = list(ATE = wls_ate,
               ATT = wls_att),
    AIPW = list(ATE = dr_ate,
                ATT = dr_att)
  )
})

cat(paste0("time: ", round(time[[3]], 2)," sec\n"))
W = res$weights

# weight check ------------------------------------------------------------


cat("\nweight sum per group\n")
print(sapply(W, function(w) tapply(w, A_, sum)))

cat("\nnbr negative weight per group\n")
print(sapply(W, function(w) tapply(w,A_, function(w) sum(w<0))))

cat("\nsum of weight above percentile per group:\n")
print(lapply(W, function(w) {
  simplify2array(lapply(quantile(w, c("90th"=.9,"95th"=.95,"99th"=.99)),
                        function(q) tapply(w, A_, function(w_grp) sum(w_grp[w_grp>q]))
  ))}
))

# weight distribution ATE -------------------------------------------------


p = data.frame(logW=log10(W$ATE), group=A_)[order(W$ATE), ] |>
  ggplot(aes(x=1:n, y=logW, colour=group)) +
  geom_point(show.legend = F) +
  ggtitle("ATE weight per group") +
  xlab("") +
  facet_wrap(vars(group), 1) +
  theme_minimal()
print(p)

p = data.frame(W=W$ATE, group=A_) |>
  ggplot(aes(W, colour = group)) +
  geom_density(show.legend = T) +
  ggtitle("ATE weight distribution per group") +
  theme_minimal()
print(p)


# weight distribution ATT -------------------------------------------------


p = data.frame(logW=log10(W$ATT), group=A_)[order(W$ATT), ] |>
  ggplot(aes(x=1:n, y=logW, colour = group)) +
  geom_point(show.legend = F) +
  ggtitle("ATT weight per group") +
  xlab("") +
  facet_wrap(vars(group), 1) +
  theme_minimal()
print(p)

p = data.frame(W=W$ATT, group=A_) |>
  ggplot(aes(W, colour = group)) +
  geom_density(show.legend = T) +
  ggtitle("ATT weight distribution per group") +
  theme_minimal()
print(p)


# estimate ATE / ATT ------------------------------------------------------


## default weight
cat(paste0(
  "\ndefault weight\n",
  "ATE =  ", sum(W$ATE * (2 * A - 1) * Y), "\n",
  "ATT =  ", sum(W$ATT * (2 * A - 1) * Y), "\n"
))

## normalize
cat(paste0(
  "\nweight normalized:\n",
  "ATE =  ", sum(norm(W$ATE, A) * (2 * A - 1) * Y), "\n",
  "ATT =  ", sum(norm(W$ATT, A) * (2 * A - 1) * Y), "\n"
))