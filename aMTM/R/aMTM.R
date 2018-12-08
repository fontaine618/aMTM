aMTM <- function(target,N,K,
                 x0,sig0,mu0,lam0,
                 adapt,global,scale,local,
                 proposal,accrate,gamma,
                 parms,
                 beta) {
   t0<-proc.time()
   # INITIAL CHECKS
   d <- length(x0)
   if(missing(parms))parms <- list(0)
   # CALL C++
   out<-aMTMsample(target,N,d,K,
              x0,sig0,mu0,lam0,
              adapt,global,scale,local,
              proposal,accrate,gamma,
              parms,
              beta)
   # OUTPUT
   names(out) <- c('X','sel','acc','alpha','mu','Sig',
                   'lam','accrate','N','K','d','Sy')
   out$time <- (proc.time()-t0)[1]
   out
}
