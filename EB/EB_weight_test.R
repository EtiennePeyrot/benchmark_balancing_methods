rm(list = ls())
graphics.off()
# set.seed(123)


# import
library(ggplot2)
source("EB_weight.R")
source("../GenerateData.R")


# input
n = 1000 # sample size
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
  W = EB(X = X, A = A)
})
cat(paste0("time: ", round(time[[3]], 2)," sec\n"))


# weight check ------------------------------------------------------------


cat("\nweight sum per group\n")
print(sapply(W, function(w) tapply(w, A_, sum)))

cat("\nnbr negative weight per group\n")
print(sapply(W, function(w) tapply(w,A_, function(w) sum(w<0))))


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