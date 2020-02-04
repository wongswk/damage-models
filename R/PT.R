load("gp-settings.rda")
library(hypergeo)
pid <- as.numeric(commandArgs()[5])
if (is.na(pid)) pid <- 1
cat(pid,"\n")


logpost <- function(pars) { llik(pars) }

temp.all <- signif(exp(seq(log(mintemp), log(maxtemp), length = nprocs)), 3)   # T from 1 to 50, nprocs threads

theta <- matrix(NA,n.iter,length(start.theta))
propscale <- 0.01 * abs(start.theta) * sqrt(temp.all[pid])
theta[1,] <- start.theta
curllik <- logpost(theta[1,])
lliklist <- c()
lliklist[1] <- curllik

for (t in 2:n.iter) {
  
  if (t %% swapint != 0) {
    # do a regular local move
    proptheta <- theta[t-1,] + rnorm(length(start.theta), 0, propscale)
    show(proptheta)
    newlik <- logpost(proptheta)
    if (is.nan(newlik))
      newlik <- 1e9
    
    lalpha <- (curllik - newlik) / temp.all[pid]
    show(lalpha)
    
    if(log(runif(1))<lalpha){
      theta[t,] <- proptheta
      curllik <- newlik
    } else{
      theta[t,] <- theta[t-1,];
    } 
    lliklist[t] <- curllik    
  }
  if (t %% swapint == (swapint - 1)) {
    # write out the parameters
    ii <- (t+1) / swapint
    cat(theta[t,], curllik, file = paste( "fit-gp/swap.",  ii, ".", pid, sep = "" ))
    
  }
  
  if (t %% swapint == 0) {
    torder <- 1:nprocs
    ii <- t / swapint
    
    if (pid == 1) {
      
      ready <- rep(F,nprocs)
      while ( sum( ready ) < nprocs ) {
        Sys.sleep(3)
        for (j in 1:nprocs) {
          ready[j] <- file.exists( paste( "fit-gp/swap.",  ii, ".", j, sep = "" ) )
        }
        show(ready)
      }
      
      curlpost <- c()
      curbetas <- matrix(nrow = nprocs, ncol = length(start.theta)+1)
      for (j in 1:nprocs) {
        proci <- scan(paste( "fit-gp/swap.",  ii, ".", j, sep = "" ))
        curlpost[j] <- proci[length(proci)]
        curbetas[j,] <- proci #proci[-length(proci)]
      }
      
      natte <- 500
      for (j in 1:natte) {
        for (kk in 2:nprocs) {
          if (runif(1) < exp((curlpost[torder[kk]] - curlpost[torder[kk-1]]) *
                             (1/temp.all[torder[kk]] - 1/temp.all[torder[kk-1]]))) {
            # Swap accepted
            swap.t <- torder[c(kk,kk-1)]
            torder[kk] <- swap.t[2]
            torder[kk-1] <- swap.t[1]
          }
        }
      }
      
      show(torder)
      
      # Dump final swaps to file.
      write.table(curbetas[torder,], file = paste( "fit-gp/swap.",  ii, ".all", sep = "" ))      
    }
    
    
    # Wait and check if the swaps index has been created.
    swapsdone <- FALSE
    while (swapsdone == FALSE) {
      swapsdone <- file.exists( paste( "fit-gp/swap.",  ii, ".all", sep = "" ) )
      Sys.sleep(2)
    }
    getPars <- read.table (  paste( "fit-gp/swap.",  ii, ".all", sep = "" ))
    
    # Store in current iterate
    theta[t,] <- unlist(getPars[pid,-ncol(getPars)])
    show(theta[t,])
    curllik <- unlist(getPars[pid,ncol(getPars)])
    lliklist[t] <- curllik
    
    save(theta,lliklist, file=paste("fit-gp/",pid,".rda",sep=""))
  }
  
}

save(theta,lliklist, file=paste("fit-gp/",pid,".rda",sep=""))


