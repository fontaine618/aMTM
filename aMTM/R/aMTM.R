aMTM <- function(target,N,K,x0,...) {
   call <- match.call()
   #-----------------------------------
   # CATCH OPTIONAL ARGUMENTS
      opt <- list(...)
      for (name in names(opt) ) assign(name, opt[[name]])
      optionalParamNames <- c("sig0","mu0","lam0","adapt","global","scale","local",
                              "proposal","accrate","gamma","parms","weight","burnin")
      unusedParams <- setdiff(names(opt),optionalParamNames)
      if(length(unusedParams))
         stop('Unused parameters: ',paste(unusedParams,collapse = ', '),'. 
              Pass target parameters as a list in the parms argument.')
   #-----------------------------------
   # INPUT CHECKS AND PREPARATION
      # extract dimension
      d <- length(x0)
      # create list of parameters if not supplied
      if(is.null(opt$parms)) opt$parms <- list(0)
      # check that dimension and parameters match for the target
      evalTarget <- tryCatch(target(x0,parms), error = function(e) "error")
      if(evalTarget == "error") stop("Cannot evaluate target. Check that x0 has the right dimension and
                 that parms contains all required parameters.")
      # Check if sig0 is supplied. Otherwise, we initialize to identity
      if(is.null(opt$sig0)){
         opt$sig0 <- sapply(seq(K),function(k)diag(d),simplify='array')
      }else{
         if(any(dim(opt$sig0)!=c(d,d,K))) stop("sig0 must be of dimension (d,d,K).")
      }
      # Check if mu0 is supplied. Otherwise, we initialize to 0
      if(is.null(opt$mu0)){
         opt$mu0 <- matrix(0,d,K)
      }else{
         if(any(dim(opt$mu0)!=c(d,K))) stop("mu0 must be of dimension (d,K).")
      }
      # Check if lam0 is supplied. Otherwise, we initialize
      if(is.null(opt$lam0)){
         opt$lam0 <- rep((2.38)^2/d,K)
      }else{
         if(dim(opt$lam0)!=K) stop("lam0 must be of dimension K.")
         if(any(opt$lam0<=0)) stop("lam0 must contain positive values.")
      }
      # Check algorithm specifications
      if(is.null(opt$adapt)) opt$adapt <- "ASWAM"
      if(!opt$adapt %in% c("AM","ASWAM","RAM",1,2,3)) warning("adapt is not one of AM, ASWAM or RAM. No adaptation is performed.")
      if(is.character(opt$adapt)) opt$adapt<-switch(opt$adapt,"AM"=1,"ASWAM"=2,"RAM"=3,0)
      if(opt$adapt == 3) opt$lam0 <- rep(1,K)
      if(is.null(opt$global)) opt$global <- F
      if(!opt$global %in% c(F,T)) stop("global must be either TRUE(1) or FALSE(0).")
      if(is.null(opt$scale)) opt$scale <- F
      if(!opt$scale %in% c(F,T)) stop("scale must be either TRUE(1) or FALSE(0).")
      if(is.null(opt$local)) opt$local <- F
      if(!opt$local %in% c(F,T)) stop("local must be either TRUE(1) or FALSE(0).")
      # Check type of proposal
      if(is.null(opt$proposal)) opt$proposal <- "common"
      if(!opt$proposal %in% c("ind","common","QMC", "EA",0,1,2,3)) stop("proposal must be one of ind(0), common(1), QMC(2), EA(3).")
      if(is.character(opt$proposal)) opt$proposal <- switch(opt$proposal,"ind"=0,"common"=1,"QMC"=2,"EA"=3)
      # Check target acceptance rate
      if(is.null(opt$accrate)) opt$accrate <- 0.234
      if(opt$accrate <= 0 || opt$accrate >= 1) stop("accrate must be between 0 and 1.")
      # Check power for adaptation step
      if(is.null(opt$gamma)) opt$gamma <- 0.7
      if(opt$gamma < 0.5 || opt$gamma > 1) stop("gamma must be between 0.5 and 1.")
      # Check weight function
      if(is.null(opt$weight)) opt$weight <- 0
      if(!opt$weight %in% c("proportional","importance",0,-1)) stop("weight must be either proportional(0) or importance(-1).")
      if(is.character(opt$weight)) opt$weight <- switch(opt$weight,"proportional"=0,"importance"=-1)
      # Check burnin
      if(is.null(opt$burnin)) opt$burnin <- 0
      if(opt$burnin < 0 || opt$burnin >= 1) stop("burnin must be between 0 and 1.")
      Nt <- ceiling(N/(1-opt$burnin))
   #-----------------------------------
   # CALL C++
   t0<-proc.time()[1]
   out<-aMTMsample(target,Nt,d,K,x0,opt$sig0,opt$mu0,opt$lam0,opt$adapt,opt$global,
                   opt$scale,opt$local,opt$proposal,opt$accrate,opt$gamma,opt$parms,opt$weight)
   time <- proc.time()[1]-t0
   names(out) <- c('X','sel','acc','mu','Sig','lam')
   #-----------------------------------
   # OUTPUT
   X <- out$X[seq(Nt-N+1,Nt),]
   acc.rate <- mean(out$acc)
   sel.prop <- table(out$sel[seq(Nt-N+1,Nt)])/N
   list(X=X,acc.rate=acc.rate,sel.prop=sel.prop,mu=out$mu,lam=out$lam,Sig=out$Sig,sel=out$sel,time=time)
}
