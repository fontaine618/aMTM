#' Adaptive Multiple-Try Metropolis algorithm
#'
#' This function performs the Adaptive Multiple-Try Metropolis algorithm as described in Fontaine and Bedard (2019). 
#' The sampling step is performed via a MTM algorithm with \eqn{K} gaussian proposals which may be correlated (see argument \code{proposal} and details) 
#' and using either weights that are proportional to the target density or importance weights. The adaptation step is performed
#' via one of AM, ASWAM or RAM updates of the selected proposal density; a global component may also be adapted at each iteration 
#' (see argument \code{global} and details) and the scale paramters may be adapted at each iteration (see argument \code{scale} and details). The AM and 
#' ASWAM update may be done using local steps rather than global steps (see argument \code{local} and details).
#' 
#' @useDynLib aMTM
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @param target Target log-density which must be vectorized, i.e. take input of dimension \eqn{K\times d}{K*d}. 
#' Additional parameters must be passed as a list trough the \code{parms} argument.
#' @param N Size of the Monte Carlo sample.
#' @param K Number of proposals in the MTM sampling.
#' @param x0 A vector of dimension \eqn{d} corresponding to the initial state of the chain.
#' @param sig0 An array of dimension \eqn{d\times d\times K}{d*d*K} containing the \eqn{K} 
#' initial covariance of the instrumental gaussian distributions. Default is \eqn{K} identity matrices.
#' @param mu0 A matrix of dimension \eqn{d\times K}{d*K} containing the \eqn{K} initial mean parameters 
#' for AM and ASWAM updates. Default is \eqn{K} zero vectors.
#' @param lam0 A vector of dimension \eqn{K} containing the \eqn{K} scale parameters. Default is 
#' \eqn{(2.38)^2/d} for AM and AWSAM updates and 1 for RAM updates.
#' @param adapt Determines the type of update of the covaraince done in the adaptation step : \dQuote{AM} or 1 performs AM updates, 
#' \dQuote{ASWAM} or 2 performs ASWAM updates, \dQuote{RAM} or 3 performs RAM updates, any other value produces no update. Default is \code{"ASWAM"}.
#' @param global Boolean value enabling the use of a global component that is updated at every iteration. Default is \code{FALSE}.
#' @param scale Boolean value enabling the adaptation of the scale parameter of the proposals that were not selected. Default is \code{FALSE}.
#' @param local Boolean value enabling the use of local update steps in AM or ASWAM updates. Default is \code{FALSE}.
#' @param proposal Determines the type of proposals used in the MTM sampling : \dQuote{ind} or 0 produces independant candidates, 
#' \dQuote{common} or 1 produces candidates from a common random vectors, \dQuote{QMC} or 2 produces candidates by a 
#' Randomized Quasi Monte Carlo procedure using a Koborov rule, \dQuote{EA} or 3 produces extremely antithetic candidates. Default is \code{"ind"}.
#' @param accrate Target acceptance rate for ASWAM and RAM updates. Default is 0.234.
#' @param gamma Power used to produce the adaptation step size : at iteration \code{n}, the stepsize is given by \code{n^-gamma}. Default is 0.7 and 
#' the range of suggested values is \eqn{(0.5,1]} to meet theoritical guarantees while values in \eqn{(0,0.5]} are accepted but may yield weird behaviour.
#' @param parms A list of paramters passed to the \code{target} function for evaluation.
#' @param weight Determines the type of weights used in the MTM sampling : \dQuote{proportional} or 0 produces weights that are proportional to the 
#' target density, \dQuote{importance} or -1 produces importance weights. Default is \code{proportional}.
#' @param burnin A value in \eqn{[0,1)} indicating the amount of burnin to be done. The Monte Carlo sample size is adjusted to yield a sample of size 
#' N after the burn-in. Default is 0.
#'
#' @return A list containing the following elements:
#' 
#' \item{\code{X}}{The MCMC sample in a matrix of dimension \eqn{N\times d}{N*d} and is of class \code{"mcmc"} to insure compatibility with
#' the \code{coda} package.}
#' \item{\code{acc.rate}}{The acceptance rate of the sample.}
#' \item{\code{sel.prop}}{The selection proportion of each proposal.}
#' \item{\code{mu}}{The final value of the mean parameters (only meaningful for AM and ASWAM updates without local steps).}
#' \item{\code{lam}}{The final value of the scale parameters.}
#' \item{\code{Sig}}{The final value of the covariance paramters.}
#' \item{\code{sel}}{A vector of length \code{N} containing the index of the selected proposal at each iteration.} 
#' \item{\code{time}}{The computing time.}
#' 
#' @details 
#' The MTM sampling with correlated candidates is described in Craiu and Lemieux (2007). The proposals all have Gaussian marginal 
#' distributions corresponding to a random walk (i.e. centered at the current state). The covariance is given by 
#' \eqn{\lambda^{(k)}\Sigma^{(k)}}{lambda^(k)Sigma^(k)} where \eqn{\lambda^{(k)}}{lambda^(k)} is a positive scale parameter and 
#' \eqn{\Sigma^{(k)}}{Sigma^(k)} is a symmetric positive-definite matrix. 
#' 
#' The correlation between the proposals are defined by one of the four following procedures. First, the candidates may be independant 
#' (\code{proposal="ind"} or \code{proposal=0}) so that the joint density of the proposal is the product of the marginals. Second, the 
#' candidates may be construted from a common random vector (\code{proposal="common"} or \code{proposal=1}) in which case a 
#' \eqn{d}-dimensional standard noraml vector is generated and the candidates are all calculated from this random 
#' vector via the natural transformation. Third, the candidates may be generated via a randomized Qausi-Monte Carlo procedure using a 
#' Koborov rule (\code{proposal="QMC"} or \code{proposal=2}) meaning that they are regularly spaced, see Fontaine and Bedard (2019) for 
#' details. Fourth, the candidates may be extremely antithetic (\code{proposal="EA"} or \code{proposal=3}) which means that the 
#' correlation is chosen to maximize the Euclidian distance between the candidates, see Fontaine and Bedard (2019) for details.
#' 
#' The weight function used in the selection and acceptation steps may be one of the following two choices. First, the weight can be chosen
#' to be proportional to the target density (\code{weight="proportional"} or \code{weight=0}), i.e. \eqn{w^{(k)}(y|x)=\pi(y)}{w^(k)(y|x)=pi(y)}.
#' Second, the weights may take the form of importance weights (\code{weight="imporatnce"} or \code{weight=-1}), i.e. 
#' \deqn{w^{(k)}(y|x)=\frac{\pi(y)}{Q^{(k)}(y|x)}.}{w^(k)(y|x)=pi(y)/Q^(k)(y|x).}
#' 
#' The adaptation of the parameters has different options. Unless no adaptation is performed, the update of the parameters of selected
#' proposal density is performed at every iteration. If \code{global=TRUE}, then the first proposal density is updated at every iteration.
#' If \code{scale=TRUE}, then all the scale parameters are adapted at every iterations. Now, the update of the parameters of a proposal may take one
#' of three forms. First, the covariance \eqn{\Sigma^{(k)}}{Sigma^(k)} can be adapted by the AM update of Haario et al. (2001) while the scale
#' parameter remains constant (\code{adapt="AM"} or \code{adapt=1}). Second, the scale parameter may be adapted to attain a target acceptance rate
#' as described by the ASWAM algorithm of Andrieu and Thoms (2008) where the covariance is adapted as in the AM update and the scale parameter is 
#' updated in the log scale to close in the gap between the mean accaptance rate and the target rate (\code{adapt="ASWAM"} or \code{adapt=2}). Third,
#' the covariance may be updated by the RAM update of Vihola et al. (2012) where a target acceptance rate is used to adapt the covariance directly
#' without updating the scale parameter (\code{adapt="RAM"} or \code{adapt=3}). In the AM and ASWAM update cases, it is possible to used local
#' steps (\code{local=TRUE}) which may help the algorithm to adjust more precisely to local strutures of the target distribution. In particular,
#' when the local moves are not used, all the covariances will often converge to the overall covariance of the target distribution.
#' 
#' In the context of adaptive MCMC using stochatic approximations to update the parameters, the adaptation stepsize \eqn{\gamma_n}{gamma_n} must 
#' satisfy some conditions. The algorithm uses \eqn{\gamma_n = n^{-\gamma}}{gamma_n=n^-gamma} for some parameter \eqn{\gamma}{gamma}. To insure 
#' good behaviour of the algorithm, we require \eqn{0<\gamma\leq 1}{0<gamma<=1}, but \eqn{0.5<\gamma\leq 1}{0.5<gamma<=1} furthur garantees the 
#' convergence of the parameters so that the transition stabilizes in the long run.
#'
#' @author Simon Fontaine, \email{fontaines@@dms.umontreal.ca}
#' 
#' @references 
#' Andrieu, C. and Thoms, J. (2008). "A tutorial on adaptive MCMC". Statistics and computing, 18:4, 343-373.
#' 
#' Craiu, R.V. and Lemieux., C. (2007). "Acceleration of the multiple-try Metropolis algorithm using
#' antithetic and stratifed sampling". Statistics and computing, 17:2, 109.
#' 
#' Fontaine, S. and Bedard, M. (2019). "An Adaptive Multiple-Try Metropolis algorithm". To be submitted.
#' 
#' Haario, H., Saksman, E., Tamminen, J. et al. (2001). "An adaptive Metropolis algorithm". Bernoulli, 7:2, 223-242.
#' 
#' Liu, J.S., Liang, F. and Wong, W.H. (2000). "The Multiple-Try Method and Local Optimization in Metropolis Sampling". 
#' Journal of the American Statistical Association, 95:449, 121-134.
#' 
#' Vihola, M. (2012). "Robust adaptive Metropolis algorithm with coerced acceptance rate". Statistics and Computing, 22:5, 997-1008.
#' 
#'
#' @examples
#' 
#' TBD
#'
#' @export
#'


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
      evalTarget <- tryCatch(target(matrix(x0,K,d,T),parms), error = function(e) "error")
      if(evalTarget == "error") stop("Cannot evaluate target. Check that x0 has the right dimension,
                 that parms contains all required parameters and that target is vectorized.")
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
      if(opt$gamma < 0 || opt$gamma > 1) stop("gamma must be between 0 and 1.")
      if(opt$gamma < 0.5) warning("We suggest using gamma between 0.5 and 1 to meet theoritical guarantees.")
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
   X <- coda::mcmc(out$X[seq(Nt-N+1,Nt),])
   acc.rate <- mean(out$acc)
   sel.prop <- table(out$sel[seq(Nt-N+1,Nt)])/N
   names(sel.prop) <- seq(K)
   list(X=X,acc.rate=acc.rate,sel.prop=sel.prop,mu=out$mu,lam=out$lam,Sig=out$Sig,sel=out$sel+1,time=time)
}
