#Function evaluating the performance of an MCMC sample
#
#ARGUMENTS
#
#X       the mcmc sample
#mlda    a lda object or a list of lda objects
#Sigma   the true covariance of the target
#f       a (multivariate) function
#f.true  the expected value of f
#
#OUTPUT
#
#msjd    the msjd of the sample
#msejd   the msejd of the sample using the true covariance
#act     the act of the sample using the true covariance
#ess     the mESS of the sample using the true covariance
#essCpu  the ESS/CPU measure
#essNbEval the ESS/NbEval measure
#dist.TV the estimated Tv distances using the mlda object(s)
#rmse    the RMSE of the estimator of E(f)
#mae     the MAE of the estimator of E(f)
#

performance.aMTM <- function(X, Sigma, f, f.true, time, NbEval, mlda){
   X <- as.matrix(X)
   #basic stats
   stats <- stats.aMTM(X, Sigma)
   ess <- stats[4]
   acc.rate <- mean(apply(diff(X),1,function(d) max(abs(d))>1e-10))
   stats <- c(stats, acc.rate=acc.rate, 
              essCpu = unname(ess)/unname(time), 
              essNbEval = unname(ess)/NbEval)
   #lda predicted regions and TV estimated distance
   if(!missing(mlda)){
      if(class(mlda)=='lda') mlda <- list(mlda)
      dist.TV <- sapply(mlda, function(obj){
         suppressWarnings(pred <- predict(obj, X))
         prop <- table(pred$class)/nrow(X)
         max(abs(prop - obj$counts/obj$N))
      })
   }else{dist.TV=NA}
   #function
   Xf <- t(apply(unname(X),1,f))
   Xf <- t(apply(Xf,2,mean))
   errorAbs <- abs(Xf-f.true)
   errorRel <- errorAbs/abs(f.true)
   #output
   list(stats=stats, errorAbs=errorAbs, errorRel=errorRel, dist.TV=dist.TV)
}

# WRAPPER FOR expand.grid TO DATA FRAME with default values for non-specifed things
experiments_expand <- function(K=1, adapt=2, global=F, scale=F, 
                        local=F, proposal=0, accrate=0.5, 
                        gamma=0.7, weight=0, burnin=0, N=1e4)  {
   data.frame(expand.grid(as.list(environment())))
}

# WRAPPER FOR PARALLEL EXPERIMENTS
#
# ARGUMENTS
#
# target       Target density, must be compatible with aMTM
# parms        arguments to be passed to target
# B            MC Replications
# x0,          initial state, required
# experiment   data.frame containing the experiments
# f            a (multivariate) function on which to evaluate the sample
# f.true       the true value of E(f)
# Sigma        The true covariance of the target
# mlda         a lda object to be used by performance (optional)
# path         a path and file name to save the results
#
# OUTPUT
#

aMTM_experiments_results <- function(target, parms, x0, B, experiments, f, f.true, Sigma, mlda, path, sig0scale){
   if(missing(x0)) stop("x0 required")
   require(doParallel) #also load foreach and parallel
   cl <- makeCluster(detectCores()-1)
   registerDoParallel(cl)
   results <- foreach(b = seq(B), 
                      .combine = list,
                      .multicombine = TRUE,
                      .verbose = F,
                      .export = c('seqLog','performance.aMTM' ),
                      .packages = c("aMTM","base","MASS"))  %dopar%  {
                         set.seed(b)
                    t(apply(experiments, 1, function(exp){
                       #produce mcmc chain
                       args <- as.list(exp)
                       args$target <- target
                       args$x0<-x0
                       args$parms <- parms
                       #initialize
                       K <- args$K
                       d <- length(x0)
                       scales <- seqLog(sig0scale, K)
                       args$mu0 <- array(0, dim = c(d,K))
                       args$lam0 <- array(2.38^2/d, dim = c(K))
                       sig0 = array(0, dim = c(d,d,K))
                       for(k in seq(K)) sig0[,,k] <- diag(d) * scales[k]
                       args$sig0<-sig0
                       #call
                       mcmc <- do.call(aMTM, args)
                       perf <- performance.aMTM(mcmc$X, Sigma, f, f.true, mcmc$time, 2*K-1, mlda=mlda)
                       unlist(perf)
                    }))
                 }
   stopImplicitCluster()
   stopCluster(cl)
   #results now contains a list of matrices, we move to array
   dim <- dim(results[[1]])
   tmpNamesStats <- dimnames(results[[1]])[[2]]
   results_array <- array(NA, dim=c(dim,B))
   for(b in seq(B)) results_array[,,b] <- results[[b]]
   #compute statistics
   results_stats <- apply(results_array, 1:2, function(x){
      c(mean = mean(x, na.rm=T), sd = sd(x, na.rm=T), min = min(x, na.rm=T), max = max(x, na.rm=T), 
        CI90 = stats::setNames(quantile(x, c(0.05, 0.95), na.rm=T, names= F),c('lower','upper')))
   })
   tmpNames <- dimnames(results_stats)[[1]]
   results_stats <- apply(results_stats, c(3,1), function(x)x)
   dimnames(results_stats) <- list(seq(nrow(experiments)), 
                                   tmpNamesStats,
                                   tmpNames
                                   )
   #first check if file already exist and update accordingly
   path <- check_file(path)
   #save results
   save(target, parms, x0, B, experiments, f, f.true, Sigma, mlda, results_stats,
        file = path)
   return(results_stats)
}

# FUCNTION TO CHECK THE FILE NAME AND UPDATES IT IF ALREADY EXISTS
check_file<-function(path, extension = ".Rdata"){
   #remove .Rdata extension
   if(substr(path, nchar(path)-5, nchar(path)) == extension) path <- substr(path, 1, nchar(path)-6)
   #check if file already exist
   if(file.exists(paste(path,extension,sep=''))){
      #check if already started counting
      if(substr(path, nchar(path)-4, nchar(path)-4) == '_'){
         n <- as.integer(substr(path, nchar(path)-3, nchar(path))) +1
      }else{
         n <- 1
      }
      path <- paste(substr(path, 1, nchar(path)-ifelse(n==1,0,4)),
                    formatC(n,width = 4, format = "d", flag = "0"),
                    sep=ifelse(n==1,'_',''))
   }
   while(file.exists(paste(path,extension,sep=''))){
      n <- as.integer(substr(path, nchar(path)-3, nchar(path))) +1
      path <- paste(substr(path, 1, nchar(path)-ifelse(n==1,0,4)),
                    formatC(n,width = 4, format = "d", flag = "0"),
                    sep=ifelse(n==1,'_',''))
   }
   #add extension
   path <- paste(path,extension,sep='')
   return(path)
}

# FUNCTION TO PRODUCE LOG SCALED SEQUENCE
seqLog <- function(lim, K){
   if(K>1){
      if(any(lim<=0))stop('lim must be positive')
      inc <- (lim[2]/lim[1])^(1/(K-1))
      seqLog <- lim[1]*inc^seq(0,K-1)
   }else{
      seqLog <- sqrt(lim[2]*lim[1])
   }
   seqLog
}
