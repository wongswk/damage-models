## Fit Gamma process model 

## Load the model file we're using
source("R/gamma-2breaks.R")  # GP with two breaks

## Settings for parallel tempering 
mintemp <- 1
maxtemp <- 10   # temperature of highest tempered chain 
nprocs <- 8    # number of parallel threads (CPU cores) to use
n.iter <- 10000 # MCMC iterations to run
swapint <- 100  # iterations between swaps

dir.create("fit-gp/")
save.image("gp-settings.rda")
## Start tempering
for (i in 1:nprocs) {
  cmd <- paste0("R --vanilla --no-save --args ", i, " < R/PT.R > fit-gp/log", i, ".txt &")
  show(cmd)
  system(cmd)
}
## MCMC results saved in "fit-gp/1.rda"
