library(deSolve)

# Number of stochastic load profiles to generate
nrep <- 100000

# Load parameters
gamma = 0.25; R_0 = 3000; alpha_d = 1.25; alpha_l = 1.5

## tau(t):  this will be stochastically generated for live load applications
zip <- function(a, b){
  idx <- order(c(seq_along(a), seq_along(b)))
  unlist(c(a,b))[idx]
}



D_d <- list()
T_s <- list()
load_s <- list()
D_s <- list()
T_e <- list()
T_p <- list()
load_p <- list()
D_e <- list()
N <- 100

for (i in 1:nrep) {

  D_d[[i]] <- rnorm(1, 1.05, 0.1)
  T_s[[i]] <- rexp(N, 0.1)
  load_s[[i]] <- rgamma(N+1, shape = 3.122, scale = 0.0481)
  D_s[[i]] <- stepfun(cumsum(T_s[[i]]), load_s[[i]])
  T_e[[i]] <- rexp(N, 1)
  T_p[[i]] <- rexp(N, 1/0.03835)
  load_p[[i]] <- rgamma(N, shape = 0.826, scale = 0.1023)
  D_e[[i]] <- stepfun(cumsum(zip(T_e[[i]], T_p[[i]])), zip(rep(0, N+1), load_p[[i]]))
  
}

save(D_d, D_s, D_e, gamma, R_0, alpha_d, alpha_l, file="rel_loads.rda")
