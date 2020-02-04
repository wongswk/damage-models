source("R/gamma-common.R")

eta <- function(t, theta, index) {
  g(sc, log_b, s1, log_s0, log_c, log_s2) %=% theta
  # Note: s1 is u, s0 is tau* 
  b <- exp(log_b)
  s0 <- exp(log_s0)
  c <- exp(log_c)
  s2 <- exp(log_s2)

  tfunc <- function(t) ifelse(t<s2, t^b, (s2)^(b-c) * t^c)
  dtfunc <- function(t) ifelse(t<s2, b*t^(b-1), (s2)^(b-c) * c * t^(c-1))
  
  
  if (T1[index] > 0 & t > T0[index] & t < T1[index] + T0[index]) {  
    ret <- sum(ifelse ( t - timevec[[index]] > 0 & levvec - s0 > 0,  tfunc(t-timevec[[index]]) * ( s1 * pmin(levvec - s0, delta_tau) ), 0))
    attr(ret,"dt") <- sum(ifelse ( t - timevec[[index]] > 0 & levvec - s0 > 0, dtfunc(t-timevec[[index]])  * ( s1 * pmin(levvec - s0, delta_tau) ), 0))
    return(ret)
  } else {
    L <- ifelse(t > T1[index], floor( (t-T1[index])*k_app[index] / delta_tau), floor( t*k_app[index] / delta_tau))
    #show(L)
    t_l <- timevec[[index]][L]
    t_u <- timevec[[index]][L+1]
    
    eta_l <- sum(ifelse ( t_l - timevec[[index]] > 0 & levvec - s0 > 0, tfunc(t_l - timevec[[index]]) * ( s1 * pmin(levvec - s0, delta_tau) ) , 0))
    eta_u <- sum(ifelse ( t_u - timevec[[index]] > 0 & levvec - s0 > 0, tfunc(t_u - timevec[[index]]) * ( s1 * pmin(levvec - s0, delta_tau) ) , 0))
    
    ret <- (eta_l * (t_u - t) + eta_u * (t - t_l) ) / (t_u - t_l)
    attr(ret,"dt") <- (eta_u - eta_l)/ (t_u - t_l)
    return(ret)
    
  }
}


llik <- function(theta) {
  
  g(sc, log_b, s1, log_s0, log_c, log_s2) %=% theta
  
  ### Constraints on breakpoints and powers
  
  if (log_b> log_c) return(1e9)
  #if (log_c> log_d) return(1e9)    
  #if (log_s3 < log_s2) return (1e9)
  #if (s0 < 0) return(1e9)
  
  llvec <- c()
  v <- 0
  for (index in 1:numData) {
    for (i in 1:length(obs[[index]])) {
      lc <- lik(obs[[index]][i],theta,index)  
      if (is.nan(lc)) { return(1e9)}
      if (lc <= 0) { return(1e9)}
      v <- v + log(lc)
      #show(c(obs[[index]][i], lc)) 
      #llvec <- c(llvec, log(lc))
      #if (lc > 100)
      #  show(c(index, i,  obs[[index]][i], lc)) 
    }
    #show(v)
  }
  
  -v  # neg llik
  #llvec
}


###  Guess of starting parameters
start.theta <- c( 0.285129704, -3.785473661,  0.000670881 , 6.461450733, -2.359606569 , 7.756773264)

### Initial crude optimization
#z1 <- optim(start.theta,llik,control=list(trace=3,maxit=10000))

