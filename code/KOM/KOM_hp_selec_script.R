rm(list = ls())
graphics.off()
# set.seed(123)

# import
source("KOM_hp_selec.R")
source("../GenerateData.R")

task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ntasks  = as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
outdir = as.character(Sys.getenv("OUTDIR"))
N_rep = as.integer(Sys.getenv("N_REPEAT"))


idx = which(rep(1:ntasks, length.out = N_rep) == task_id)
N = length(idx) # number of iteration for current task

dimname = list(n = c(250, 500, 1000, 2000),
               prop_trt = c("low", "moderate", "high"),
               confdg_lvl = c("low", "moderate", "high"))

# create output (now storing lambda, scale_feat, scale_arm)
hp = array(
  data = NA_real_,
  dim  = c(sapply(dimname, length), 2, 2, N),
  dimnames = c(dimname, grp=list(c("ctrl", "trt")), hp=list(c("lambda", "scale"))
  )
)

pars = do.call(expand.grid, c(dimname, stringsAsFactors = FALSE))

for (row in 1:nrow(pars)) {
  par = unlist(pars[row, ])
  
  for (i in seq_along(idx)) {
    set.seed(idx[i])
    data = GenerateData(
      n = as.integer(par[1]),
      prop_trt  = par[2],
      confdg_lvl  = par[3]
    )[1:3]
    
    res = do.call(KOM_hp, data)
    
    hp[par[1], par[2], par[3], "ctrl", , i] = res$ctrl
    hp[par[1], par[2], par[3], "trt" , , i] = res$trt
  }
}


save(hp, file=paste0(outdir,"/res/",task_id,".RData"))