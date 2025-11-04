rm(list = ls())
graphics.off()
set.seed(125)

source("KOM_hp_selec.R")
source("../GenerateData.R")


# input
n = 2000 # sample size
prop_trt = "moderate" # c("low", "moderate", "high")
confdg_lvl = "moderate" # c("low", "moderate", "high")


# generate data
data = GenerateData(n = n, prop_trt = prop_trt, confdg_lvl = confdg_lvl)
X = as.matrix(data$X)
A = as.logical(data$A)
A_ = ifelse(A,"treated","control")
Y = data$Y

# compute weights
time = system.time({
  hp = KOM_hp(X = X, A = A, Y = Y)
})

cat(paste0("\ntime = ", round(time[[3]], 2), " sec\n"))
cat(paste0("\nbest hyper parameter for n = ", n, "\n"))
print(simplify2array(hp))