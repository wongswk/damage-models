## All gamma process models source this file (Jan. 3, 2020)

library(hypergeo)

source("R/loadData.R")
source("R/bunches.R")
numData <- length(obs)

delta_tau <- 20
levvec <- seq(delta_tau,tau_max,by=delta_tau)  

timevec <- list()
for (i in 1:numData) {
  if (T1[i] == 0)
    timevec[[i]] <- levvec / k_app[i]   # time when load began to exceed that level
  else
    timevec[[i]] <- ifelse(levvec <= tauc[i], levvec / k_app[i], T1[i] + levvec / k_app[i])  
}

lik <- function(t,theta,index) {
  
  foo <- eta(t,theta,index)
  eta_t <- as.numeric(foo)
  eta_dt <- attr(foo,"dt")
  
  #show(c(index,t,eta_t,eta_dt))
  if (is.nan(eta_t)) { return(0) }
  if (eta_t > 170 || eta_t ==0) { return(0) }
  
  sc <- theta[1]
  eta_dt * (digamma(eta_t) - log(1/sc)) * pgamma(1/sc,eta_t)   +  eta_dt /eta_t^2 /gamma(eta_t) * (1/sc)^eta_t * genhypergeo(c(eta_t,eta_t),c(eta_t+1,eta_t+1),-1/sc)
}

tCDF <- function(t,theta,index) {
  eta_t <- eta(t,theta,index)
  if (is.na(eta_t))
    eta_t <- 1e9
  sc <- theta[1]
  1-pgamma(1/sc,eta_t)
  
}
