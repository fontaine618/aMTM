#' @name stats.aMTM
#' @export stats.aMTM
#' @method stats aMTM
#' @title Statistics of a MCMC chain
#'
#' @description Summary statistics for the output chain of an aMTM object.
#' 
#' @param X A matrix coresponding to the outpout of a MCMC algorithm.
#' 
#' @param cov A covariance matrix to compute the MSEJD with. Default is \code{NULL} and uses the sample covariance.
#' 
#' @return A vector containing the following statistics:
#' 
#' \item{\code{msejd}}{The Mean Sqaured Euclidian Jumping Distance.}
#' \item{\code{msjd}}{The Mean Sqaured Jumping Distance (using the sample variance in the Mahalanobis distance).}
#' \item{\code{act}}{The Frobenius norm of the multivariate ACT of the chain.}
#' \item{\code{ess}}{The multivariate ESS of the chain as described by Vats et al. (2015).}
#' 
#' 
#' @author Simon Fontaine, \email{simfont@@umich.edu}
#' 
#' @seealso \link{aMTM}.
#'
#' @references 
#' 
#' Fontaine, S. and Bedard, M. (2019). "An Adaptive Multiple-Try Metropolis algorithm". To be submitted.
#' 
#' Vats, D., Flegal, J. M., and, Jones, G. L. (2015). "Multivariate Output Analysis for Markov chain Monte Carlo". arXiv preprint arXiv:1512.07713.
#' 
#' 
#' @examples
#' \dontrun{
#' library(aMTM)
#' # Banana log-density with parameter B and a
#' p <- function(x, p) apply(x,1,function(x) -x[1]^2/(2*p$a^2) - 1/2*(x[2]+p$B*x[1]^2-p$B*p$a^2)^2)
#' # setup
#' set.seed(1)
#' N<-1e5;K<-3
#' B<-0.04;a<-8
#' # aMTM sampling with ASWAM update
#' mcmc <- aMTM(target=p, N=N, K=K, x0=c(0,0), parms=list(a=a,B=B), burnin=0.1)
#' 
#' stats.aMTM(mcmc$X)
#' 
#' stats.aMTM(mcmc$X, diag(1:2))
#'}
#'

stats.aMTM <- function(X, cov = NULL){
   #mean squared (euclidian) jump distance
   dif <- diff(X)
   msjd <- mean(sqrt(apply(dif^2,1,sum)))
   if(is.null(cov)) {
      Sigma <- var(X)
   }else{
      if(any(dim(cov) != ncol(X))) stop("cov must have the same dimension as X")
      if(any(abs(cov - t(cov))>1e-10)) stop("cov must be symmetric")
      if(any(eigen(cov)$values<1e-10)) stop("cov must be positive definite")
      Sigma <- cov
   }
   Sigmainv <- solve(Sigma)
   msejd <-  mean(sqrt( apply(dif, 1, function(row) row %*% Sigmainv %*% row) ))
   #autocorrelation time of the mean
   SigmaP <- mcmcse::mcse.multi(X)$cov
   S <- t(chol(Sigma))
   Sinv <- solve(S)
   ACT <- Sinv %*% SigmaP %*% t(Sinv)
   act <- sqrt(sum(ACT^2))#frobenius of ACT
   #multivariate ESS
   # ess <- mcmcse::multiESS(X, cov) does not allow user supplied covariance...
   p <- ncol(X)
   N <- nrow(X)
   det.var.p <- prod(eigen(Sigma, only.values = TRUE)$values^(1/p))
   det.covmat.p <- prod(eigen(SigmaP, only.values = TRUE)$values^(1/p))
   ess <- N * (det.var.p/det.covmat.p)
   #output
   round(c(msejd = msejd, msjd = msjd,act = act, ess=ess),3)
}


