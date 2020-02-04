## Gamma process reliability assessments

### Load model file used and MCMC samples
source("R/gamma-2breaks.R")
load("fit-gp/1.rda")

# Define phi 
phi <- 1.2
#phi <- as.numeric(commandArgs()[5])
show( paste0("phi is ", phi))

# Load pre-saved profiles
load("rel_loads.rda")
nrep <- length(D_d)

# Pick a thinned subset of thetas to use for reliability calculations
theta_thin <- theta[seq(nrow(theta)/2+1,nrow(theta), length.out = 500),]

failmat <- matrix(NA, nrow = nrow(theta_thin), ncol = nrep)

for (ii in 1:nrep) {
  show(ii)

  tau <- function(t) {
    # express t in years
    phi*R_0*(gamma*D_d[[ii]] + D_e[[ii]](t) + D_s[[ii]](t))/(gamma*alpha_d + alpha_l)
  }
  
  delta <- 0.001
  times <- seq(0, 50, delta)  # times for which to solve alpha(t).
  #plot(times, tau(times), type = "l", xlab = "Year", ylab = "Load(psi)")
  tau_t <- tau(times)
  
  nlevels <- floor(max(tau_t)/delta_tau)
  
  load_prof <- matrix(nrow=nlevels, ncol=length(times))
  load_prof[,1] <- 0
  
  for (i in 1:nlevels) {
    #show(i)
    for (j in 2:length(times)) {
      load_prof[i,j] <- load_prof[i,j-1] + delta * (tau_t[j] > delta_tau * i)
    }
  }
  
  eta_loadex <- function(t, theta) {  ## currently set up for 2 breakpoints
    g(sc, log_b, s1, log_s0, log_c, log_d, log_s2, log_s3) %=% theta
    # Note: s1 is u, s0 is tau* 
    b <- exp(log_b)
    s0 <- exp(log_s0)
    c <- exp(log_c)
    d <- exp(log_d)
    s2 <- exp(log_s2)
    s3 <- exp(log_s3)
    
    t_ind <- round(t/delta,0)
    tv <- load_prof[,t_ind] * 8760
    
    
    tfunc <- function(t) ifelse(t<s2, t^b, ifelse(t<s3,(s2)^(b-c) * t^c, (s2)^(b-c) * s3^(c-d) * t^d))
    dtfunc <- function(t) ifelse(t<s2, b*t^(b-1), ifelse(t<s3,(s2)^(b-c) * c * t^(c-1), (s2)^(b-c) * s3^(c-d) * d*t^(d-1)))
  
    ret <- sum(ifelse ( tv > 0 & levvec[1:nlevels] - s0 > 0,  tfunc(tv) * ( s1 * pmin(levvec[1:nlevels] - s0, delta_tau) ), 0))
    attr(ret,"dt") <- sum(ifelse ( tv > 0 & levvec[1:nlevels] - s0 > 0, dtfunc(tv)  * ( s1 * pmin(levvec[1:nlevels] - s0, delta_tau) ), 0))
    
    ret  
  }
  
  eta_loadex_all <- c()

  for (i in 1:nrow(theta_thin)) {
    eta_loadex_all[i] <- eta_loadex(50,theta_thin[i,])
  }
  
  fail50 <- 1 - pgamma ( 1 / theta_thin[,1], eta_loadex_all)

  failmat[,ii] <- fail50
}

save(failmat, file = paste0("GPprob_", phi, ".rda"))
