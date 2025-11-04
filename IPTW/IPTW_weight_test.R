rm(list = ls())
graphics.off()
# set.seed(123)


# import
library(ggplot2)
source("IPTW_weight.R")
source("confidence_interval.R")
source("../GenerateData.R")


# function
norm = function(w) w / ifelse(A, sum(w[A]), sum(w[!A]))
trunc = function(w,p) pmax(w,quantile(w,p))


# input
n = 2000 # sample size
prop_trt = "moderate" # c("low", "moderate", "high")
confdg_lvl = "high" # c("low", "moderate", "high")


# generate data
data = GenerateData(n = n, prop_trt = prop_trt, confdg_lvl = confdg_lvl)
X = data$X
A = as.logical(data$A)
A_ = ifelse(A,"treated","control")
Y = data$Y


# compute weights
time = system.time({
  ps = IPTW_logreg_weight(X = X, A = A)
  
  
  # Propensity block for variance sandwich estimator
  ps_block = propensity_blocks_logit_mle(A = A, X = X, model.fit = ps$model)
  
  wls_ate = wls_ci(Y, A, ps_block = ps_block, target = "ATE", alpha = 0.05)
  wls_att = wls_ci(Y, A, ps_block = ps_block, target = "ATT", alpha = 0.05)
  
  # Prognostic models
  outcome_block = outcome_blocks_binomial(Y = Y, A = A, X = X)
  
  # AIPW
  dr_ate = dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block, target = "ATE", alpha = 0.05)
  dr_att = dr_ci(A = A, outcome_block = outcome_block, ps_block = ps_block, target = "ATT", alpha = 0.05)
  
  estimate = list(
    WLS = list(ATE = wls_ate,
               ATT = wls_att),
    AIPW = list(ATE = dr_ate,
                ATT = dr_att)
  )
})
W = ps$W
cat(paste0("time: ", round(time[[3]], 2)," sec\n"))


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
  ggtitle("ATE weight per group") +
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
  "ATE =  ", sum(norm(W$ATE) * (2 * A - 1) * Y), "\n",
  "ATT =  ", sum(norm(W$ATT) * (2 * A - 1) * Y), "\n"
))

## truncate
p = c("90th"=.9,"95th"=.95,"99th"=.99)
for (i in seq_along(p)) {
  # without normalization
  W_trunc = list(ATE = trunc(W$ATE, p[i]), ATT = W$ATT)
  W_trunc$ATT[!A] = trunc(W$ATT[!A], p[i])
  
  cat(paste0(
    "\nweight truncation at ",names(p)[i],"-percentile:\n",
    "ATE =  ", sum(W_trunc$ATE * (2 * A - 1) * Y), "\n",
    "ATT =  ", sum(W_trunc$ATT * (2 * A - 1) * Y), "\n"
  ))
  
  # with normalization
  cat(paste0(
    "\nweight truncation at ",names(p)[i],"-percentile and normalized:\n",
    "ATE =  ", sum(norm(W_trunc$ATE) * (2 * A - 1) * Y), "\n",
    "ATT =  ", sum(norm(W_trunc$ATT) * (2 * A - 1) * Y), "\n"
  ))
}
