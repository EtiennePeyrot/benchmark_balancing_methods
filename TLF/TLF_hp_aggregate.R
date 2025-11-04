N = 50
array_len = 50
dir = paste0(choose.dir(),"\\")

task_idx = rep(1:array_len, length.out=N)

task_success = sapply(1:array_len, function(i) {
  j = data.table::first(which(task_idx==i))
  !inherits(try(load(paste0(dir,j,".RData")), silent=T), "try-error")
})

N_ = which(task_idx %in% which(task_success))
load(paste0(dir,data.table::first(N_),".RData"))

hp = array(data = NA_real_,
           dim = c(dim(res)[-8],length(N_)),
           dimnames = c(dimnames(res)[-8],list(seed=N_)))

for (i in which(task_success)) {
  idx = which(task_idx == i)
  for (j in idx) {
    load(paste0(dir,j,".RData"))
    hp[ , , , , , , , as.character(j)] = res
  }
}

save(hp, file = "hp.RData")