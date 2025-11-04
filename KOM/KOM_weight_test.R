rm(list = ls())
graphics.off()
set.seed(123)


# import
library(ggplot2)
source("KOM_weight.R")
source("confidence_interval.R")
source("../GenerateData.R")

# pre computed hyper parameters
hp = array(
  data = c(0.7786532, 0.6011859, 0.5868907, 0.6063836, 0.98312  , 0.7602185, 0.7513465, 0.6081164, 1.208586 ,
           1.193911 , 0.8165049, 0.6628186, 0.4031446, 0.4091343, 0.3962227, 0.3755389, 0.4671336, 0.50134  ,
           0.4401676, 0.4012346, 0.6484486, 0.6394642, 0.499708 , 0.4825143, 0.2631844, 0.2708328, 0.2316773,
           0.2707671, 0.3153764, 0.2790501, 0.2607347, 0.2186607, 0.2925336, 0.3830627, 0.3129639, 0.2682089,
           0.9140153, 0.6580625, 0.4573556, 0.455861 , 0.7767795, 0.547296 , 0.5210961, 0.4919411, 0.6802308,
           0.5107976, 0.5624083, 0.5249304, 0.330268 , 0.3143082, 0.275806 , 0.3108134, 0.3430622, 0.3018467,
           0.3180238, 0.3422981, 0.3159365, 0.3556322, 0.360707 , 0.4001348, 0.1803495, 0.1682394, 0.2238222,
           0.2979886, 0.1733366, 0.1846091, 0.2460524, 0.3522965, 0.213528 , 0.2231472, 0.283161 , 0.3637716,
           0.1336119, 0.1235754, 0.1349632, 0.1354374, 0.1256994, 0.1236862, 0.1412542, 0.129654 , 0.1431392,
           0.1316681, 0.1239864, 0.1149471, 0.1816869, 0.1783617, 0.1793874, 0.179703 , 0.1933758, 0.1882647,
           0.1778932, 0.1801886, 0.2087735, 0.1694415, 0.1976844, 0.1807149, 0.227526 , 0.2202755, 0.2315285,
           0.2442348, 0.2441049, 0.228016 , 0.2363547, 0.2230652, 0.2685443, 0.2527799, 0.2618525, 0.2449143,
           0.133382 , 0.1230235, 0.1056677, 0.1106035, 0.1306839, 0.1169748, 0.1324153, 0.1272273, 0.1187675,
           0.1138824, 0.1415419, 0.1293276, 0.1742895, 0.1666444, 0.1494678, 0.1512446, 0.1577851, 0.1452945,
           0.1557642, 0.1695934, 0.153698 , 0.1676022, 0.167052 , 0.1857751, 0.1732412, 0.1790957, 0.1962752,
           0.2261315, 0.1736528, 0.1865749, 0.2136447, 0.2536154, 0.1885792, 0.2001411, 0.2301916, 0.2650601),
  dim = c(4, 3, 3, 2, 2),
  dimnames = list(n = c(250, 500, 1000, 2000),
                  prop_trt = c("low", "moderate", "high"),
                  confdg_lvl = c("low", "moderate", "high"),
                  grp = c("ctrl", "trt"),
                  hp = c("lambda", "scale")))


# input
n = 250 # sample size
prop_trt = "moderate" # c("low", "moderate", "high")
confdg_lvl = "high" # c("low", "moderate", "high")


# generate data
data = GenerateData(n = n, prop_trt = prop_trt, confdg_lvl = confdg_lvl)
X = as.matrix(data$X)
A = as.logical(data$A)
A_ = ifelse(A,"treated","control")
Y = data$Y

# pre computed hyper parameters

lambda = as.list(hp[as.character(n),prop_trt,confdg_lvl,,"lambda"])
scale  = as.list(hp[as.character(n),prop_trt,confdg_lvl,,"scale"])

# compute weights
time = system.time({
  W = KOM_weight(X = X, A = A, lambda = lambda, scale = scale)
  W_ = lapply(W, pmax, 0) # prevent precision error
  
  wls_ate = wls_ci_optimist(Y, A, W = W_$ATE, alpha = 0.05)
  wls_att = wls_ci_optimist(Y, A, W = W_$ATT, alpha = 0.05)
  
  # Prognostic models
  outcome_block = outcome_blocks_binomial(Y = Y, A = A, X = X)
  
  # AIPW
  dr_ate = dr_ci_optimist(A = A, outcome_block = outcome_block, W = W_$ATE, target = "ATE", alpha = 0.05)
  dr_att = dr_ci_optimist(A = A, outcome_block = outcome_block, W = W_$ATT, target = "ATT", alpha = 0.05)
  
  estimate = list(
    WLS = list(ATE = wls_ate,
               ATT = wls_att),
    AIPW = list(ATE = dr_ate,
                ATT = dr_att)
  )
})
cat(paste0("time: ", round(time[[3]], 2)," sec\n"))

# weight check ------------------------------------------------------------


cat("\nweight sum per group\n")
print(sapply(W, function(w) tapply(w, A_, sum)))

cat("\nnbr negative weight per group\n")
print(sapply(W, function(w) tapply(w,A_, function(w) sum(w<0))))

cat("\nsum negative weight per group\n")
print(sapply(W, function(w) tapply(w,A_, function(w) sum(pmin(w,0)))))

# set negative weight to 0
if(any(sapply(W, function(w) w<0)))
  message("negative weight set to 0")
W = lapply(W, pmax, 0)

# weight distribution ATE -------------------------------------------------


p = data.frame(logW=log10(W$ATE), group=A_)[order(W$ATE), ] |>
  ggplot(aes(x=1:n, y=logW, colour=group)) +
  geom_point(show.legend = F) +
  ggtitle("ATE weight per group") +
  xlab("") +
  facet_wrap(vars(group), 1) +
  theme_minimal()
print(p)

mask = W$ATE > 1e-7
cat("\nnbr weight ATE above / below threshold\n")
aggregate(W$ATE, list(A_, mask), length, drop=F)[ ,3] |>
  matrix(2,2,F,list(c("control", "treated"), c("W < 1e-7", "W > 1e-7"))) |>
  print()

cat("\nsum weight ATE above / below threshold\n")
aggregate(W$ATE, list(A_, mask), sum, drop=F)[ ,3] |>
  matrix(2,2,F,list(c("control", "treated"), c("W < 1e-7", "W > 1e-7"))) |>
  print()

p = data.frame(W=W$ATE, group=A_)[mask, ] |>
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

mask = W$ATT > 1e-7
cat("\nnbr weight ATT above / below threshold\n")
aggregate(W$ATT, list(A_, mask), length, drop=F)[ ,3] |>
  matrix(2,2,F,list(c("control", "treated"), c("W < 1e-7", "W > 1e-7"))) |>
  print()

cat("\nsum weight ATT above / below threshold\n")
aggregate(W$ATT, list(A_, mask), sum, drop=F)[ ,3] |>
  matrix(2,2,F,list(c("control", "treated"), c("W < 1e-7", "W > 1e-7"))) |>
  print()

p = data.frame(W=W$ATT, group=A_)[mask, ] |>
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