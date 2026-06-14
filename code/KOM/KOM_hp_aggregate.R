rm(list = ls())


N = 50
ntasks = 50
dir = "kom_hp_2245954/res/"

my_load = function(i) load(paste0(dir, i, ".RData"), envir = globalenv())

idx = rep(1:ntasks, length.out = N)

my_load(1)
res = array(data = NA_real_,
            dim = c(dim(hp)[1:(length(dim(hp))-1)],N),
            dimnames = dimnames(hp))


for(task_id in 1:ntasks) {
  my_load(task_id)
  id = which(idx == task_id)
  res[,,,,,id] = hp
}
hp = res; rm(res)

save(hp, file="kom_hp.RData")