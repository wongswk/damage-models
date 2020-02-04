### Fit US Model   (Jan 2, 2020)

tau0 <- 6468 # Median short-term strength

##  Setup
source("R/loadData.R")
numData <- length(obs)

source("R/bunches.R")
EOS <- list()
for (i in 1:numData) {
  EOS[[i]] <- EnvStats::evNormOrdStats(length(obs[[i]]))
}

wNLS <- function(logT, R, k, tauc, Tp0, Tp1, case, wgt, a, b, w) {

  ewR <- exp(w * R)
  alpha_T1 <- ifelse(case ==2, (Tp1 - Tp0) * exp ( -a + b * tauc / ewR) + exp(-a) * ewR / (b*k) * ( exp( b * k * Tp0 / ewR) - 1), 0)
  res <- ifelse( case == 0, logT - log(ewR / (b * k) * log ( exp(a) * b * k / ewR + 1)),
          ifelse ( case == 1, (logT - log(tauc/k - ewR / (b*k) + exp( -b * tauc / ewR) * ( ewR / (b*k) + exp(a)))) * wgt, 
                  logT - log( ewR / (b * k) * log ( exp(a) * b * k / ewR * (1 - alpha_T1)+ 1)   ) )   #case = 2
                  
  )
  
  if(sum(alpha_T1>1)>0) message(paste("alpha_T1 exceeded in sample# ", which(alpha_T1>1)))  ## bad cases
  ifelse( is.nan(res), (alpha_T1-1)*1000, res)
  
}

alldat <- data.frame(T = (unlist(obs[1:numData])),
                     R = unlist(EOS[1:numData]),
                     k = unlist(sapply(1:numData, function(i) rep(k_app[i], length(obs[[i]])))),
                     i = unlist(sapply(1:numData, function(x) rep(x, length(obs[[x]])))),
                     tauc = unlist(sapply(1:numData, function(x) rep(tauc[x], length(obs[[x]])))),
                     Tp0 = unlist(sapply(1:numData, function(x) rep(T0[x], length(obs[[x]])))),
                     Tp1 = unlist(sapply(1:numData, function(x) rep(T1[x], length(obs[[x]])))),
                     wgt = 1
                     )

alldat$logT <- with(alldat, {ifelse( T > T1[i] & T1[i] > 0, log(T - T1[i]), log(T) )})
alldat$case <- with(alldat, {ifelse( T > T1[i] & T1[i] > 0, 2, ifelse(T > T0[i] & T0[i] > 0, 1, 0))})

fit1 <- nls( ~ wNLS(logT, R, k, tauc, Tp0, Tp1, case, wgt, a, b, w), data = alldat, start=list(a=20.665529538, b = 0.000183925, w = 0.803524486), trace=T, control=list(maxiter=100,printEval = T))
fit1s <- summary(fit1)
itcursigma <- 1e9

alldat$wgt <- ifelse(alldat$case==1, 1/fit1s$coefficients[2]/alldat$tauc , 1)
 
## Exclude any problematic cases here if necessary


### Iterate until convergence
while( abs(itcursigma - fit1s$sigma) > 1e-8) {
  itcursigma <- fit1s$sigma
  fit1 <- nls( ~ wNLS(logT, R, k, tauc, Tp0, Tp1, case, wgt, a, b, w), data = alldat, start=list(a = fit1s$coefficients[1], b = fit1s$coefficients[2], w = fit1s$coefficients[3]), trace=T, control=list(maxiter=1000,printEval = T))
  fit1s <- summary(fit1)
  alldat$wgt <- ifelse(alldat$case==1, 1/fit1s$coefficients[2]/alldat$tauc , 1)
}

ger_a <- fit1s$coefficients[1]
ger_b <- fit1s$coefficients[2]  # this is actually B' (we use it throughout in USADM as well)
ger_w <- fit1s$coefficients[3]

## Write out parameter estimate
write.table( cbind(ger_a,ger_b,ger_w), file="US_theta.csv", row.names = F, sep=',', col.names = F)
       