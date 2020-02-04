#### Load data files and set up variables

obs <- list()
obs[[1]] <- scan("data/ramp_388440.txt")
obs[[2]] <- scan("data/RCR_4500_1Y.txt")
# additional datasets can be added by extending this list

k_s <- 388440  # standard ramp-loading rate
tau_max <- 20000 # for gamma process model, the maximum strength of the lumber population

## The following vectors should match the order of datasets loaded above.

k_app <- c(388440, k_s) # ramp-loading rates applied (may have rates different from k_s for ramp load data)
T1 <- c(0, 8760)  # constant load duration (use 0 for ramp load data)
tauc <- c(Inf, 4500)  # constant load level (use Inf for ramp load data)
T0 <- c(0, 4500/388440) # time to reach constant load level (use 0 for ramp load data)
