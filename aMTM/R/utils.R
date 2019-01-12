#' Plot the chain of an aMTM object
#'
#' An object of class aMTM contains the resulting chain \code{X} as well as the selected proposal at each iteration and the resulting
#' proposal densities. This function produces various plots to inspect the output of an aMTM object.
#' 
#' @param aMTMobject An object of class \code{aMTM}.
#' @param vars A vector of indices of the variables to plot. A maximum of 10 variables is accepted and the default 
#' is all variables of \code{X} up to the first 10.
#' @param type Whether to plot points (\code{'p'}), line segments between points (\code{'l'}) or both (\code{'l'}) 
#' where the upper triangular matrix contains the trace plots and the lower diagonal matrix contains the scatter plot. 
#' Only applies when \code{pairs==TRUE}; default is \code{'p'}.
#' @param color Whether to plot the points/lines colored with the corresponding proposal density. 
#' In the case of points, the color represent the density used to get to the point; in the case of lines, 
#' the color represent the density used for the jump. Default is \code{TRUE}.
#' @param pairs Whether to plot the variables by pairs or not. Default is \code{TRUE}.
#' @param prop.density Whether to add the resulting densities to the plots or not. Default is \code{TRUE}.
#' @param ... Additionnal graphical parameters.
#'
#' @details 
#' When \code{pairs=TRUE}, the 2-D marginal plots of \code{X} for each pair of variables in \code{vars} are plotted. On the diagonal, the
#' 1-D marginal densities are plotted. When \code{pairs=FALSE}, the trace plot and the 1-D marginal plot are produced for each varaibles in 
#' \code{vars}.
#'
#' @author Simon Fontaine, \email{fontaines@@dms.umontreal.ca}
#' 
#'
#' @examples
#' 
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
#' #plot the pairs showing color
#' plot.aMTM(mcmc, color=T)
#' #plot the pairs showing jumps and color
#' plot.aMTM(mcmc, type='l', color=T)
#' #plot the marginals with colors
#' plot.aMTM(mcmc, pairs=F, color=T)
#'
#' @export
#'


plot.aMTM <- function(aMTMobject, vars, type, color, pairs, prop.density, ...){
   #type contains the type of plot to produce
   #either points or lines
   #color=T trigger the use of colors to represent the candidats
   #pairs = T triggers the plotting of the pairs rather than each trace plot
   #vars is the identifyer for which vars to plot
   #  if only lenght one : only marginal
   #  if length >1 then pairs with marginals on the diagonal
   #density plots the proposal densities at the last iteration
   #... are additionnal graphical parameters
   
   X <- aMTMobject$X
   K <- length(aMTMobject$sel.prop)
   #check that vars is subset of all possible vars
   if(missing(vars)) vars <- seq(ncol(X))
   if(any(!vars %in% seq(ncol(aMTMobject$X)))) stop('vars must be a subset of the possible range of dimensions')
   if(length(vars)>10){
      vars <- vars[1:10]
      warning('vars must be of length 10 or less. The first ten indices plotted only.')
   }
   vars <- vars[order(vars)]
   np <- length(vars)
   N <- nrow(X)
   #check other input
   if(missing(type)) type <- 'p'
   if(!type %in% c('l','p','b')) stop('type must be either l, p or b')
   if(missing(color)) color <- T
   if(!color %in% c(T,F)) stop('color must be either T or F')
   if(missing(pairs)) pairs <- T
   if(!pairs %in% c(T,F)) stop('pairs must be either T or F')
   if(missing(prop.density)) prop.density <- T
   if(!prop.density %in% c(T,F)) stop('prop.density must be either T or F')
   
   #precompute the means for density
   if(prop.density){
      #this is a p*K matrix
      mu <- sapply(seq(K), function(k){
         ids <- which(aMTMobject$sel == k)-1
         if(ids[1] == 0)ids <- ids[-1]
         colMeans(X[ids,])
      }, simplify = 'array')
   }
   #precompute colors
   if(color){
      cols <- aMTMobject$sel
   }else{
      cols <- rep(1,N)
   }
   #ranges in a 2*p matrix
   p <- ncol(X)
   ranges <- sapply(seq(p), function(j) range(X[,j]), simplify='array')
   
   #the two cases are given by pairs
   if(pairs){
      #we plot pairs, so we initialize the graph to a np x np grid
      par(mfrow=c(np,np), mar = c(0,0,0,0), oma = c(3,3,3,3), ...)
      #R fills the plots by rows
      for(i in vars){ #row i 
         for(j in vars){ #row i column j
            if(j == i){
            # on the diagonal, we plot the marginal
               plot(NA, xlim = ranges[,j], ylim = ranges[,j], 
                    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
               d <- stats::density(X[,j])
               d$y <- d$y * diff(ranges[,j])/ (2*max(d$y))
               lines(x = d$x, y=d$y + ranges[1,j])
               lines(y = d$x, x=ranges[2,j] - d$y)
               #text(x =  ranges[1,j], y =  ranges[2,j], labels = paste('Var',j), cex = 2)
               legend(x='topleft', legend = paste('Var',j), cex = 1, bty = 'n')
            }else{
            # off diagonal we plot the pairs
               plot(NA, xlim = ranges[,j], ylim = ranges[,i], 
                    xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
            # we plot the ticks on even top left, odd bottom right
               if(i==vars[1] & j%%2 == 0) axis(side=3)
               if(j==vars[1] & i%%2 == 0) axis(side=2)
               if(i==vars[np] & j%%2 == 1) axis(side=1)
               if(j==vars[np] & i%%2 == 1) axis(side=4)
               if(type=='p' | (type == 'b' & j<i) )points(X[,c(j,i)], col = cols)
               if(type=='l' | (type == 'b' & j>i) ){
                  segments(x0 = as.vector(X[-N,j]),
                           x1 = as.vector(X[-1,j]),
                           y0 = as.vector(X[-N,i]),
                           y1 = as.vector(X[-1,i]),
                           col = cols[-1])
               }
               if(prop.density){
                  #for(k in seq(K)) points(x=mu[j,k], y=mu[i,k], pch=19, col=1, cex=1.5)
                  #for(k in seq(K)) points(x=mu[j,k], y=mu[i,k], pch=19, col=k, cex=1.2)
                  for(k in seq(K)) mixtools::ellipse(mu[c(j,i),k], aMTMobject$Sig[c(j,i),c(j,i),k]*aMTMobject$lam[k],
                                                     alpha = 0.1, col=k)
               }
            }#end off diagonal
         }#end j loop
      }#end i loop
   }else{
      #pairs == FALSE
      lay <- rep(c(0,0,1),np)+floor(seq(0,3*np-1)/3)*2+1
      layout(matrix(lay, nrow = np, ncol = 3, byrow = TRUE))
      par(mar = c(3,3,3,3), oma = c(0,0,0,0),...)
      for(i in vars){
         #loop through the vars
         #first the trace plot
         plot(NA, xlim = c(1,N), ylim = ranges[,i],  xlab = '', ylab = '', main = paste('Trace of var', i))
         segments(x0 = seq(N-1),
                  x1 = seq(N-1)+1,
                  y0 = as.vector(X[-N,i]),
                  y1 = as.vector(X[-1,i]),
                  col = cols[-1])
         #second the density 
         d <- stats::density(X[,i])
         plot(NA, xlim = range(d$y), ylim = ranges[,i],  xlab = '', yaxt= 'n', ylab = '', xaxt = 'n',
              main = paste('Density of var', i))
         lines(x=d$y, y=d$x)
         if(prop.density){
            for(k in seq(K)) lines(y=d$x, 
                                   x=dnorm(d$x, mu[c(i),k], sqrt(aMTMobject$Sig[i,i,k]*aMTMobject$lam[k])),
                                   col = k, lty =3 )
         }
      }
   }
}



#' Statistics of a MCMC chain
#'
#' Summary statistics for the output chain of an aMTM object.
#' 
#' @param X A matrix coresponding to the outpout of a MCMC algorithm.
#'
#' @return A vector containing the following statistics:
#' 
#' \item{\code{msejd}}{The Mean Sqaured Euclidian Jumping Distance.}
#' \item{\code{msjd}}{The Mean Sqaured Jumping Distance (using the sample variance in the Mahalanobis distance).}
#' \item{\code{act}}{The Frobenius norm of the multivariate ACT of the chain.}
#' \item{\code{ess}}{The multivariate ESS of the chain as described by Vatset al. (2015).}
#' 
#' 
#' @author Simon Fontaine, \email{fontaines@@dms.umontreal.ca}
#' 
#'
#' @references 
#' 
#' Fontaine, S. and Bedard, M. (2019). "An Adaptive Multiple-Try Metropolis algorithm". To be submitted.
#' 
#' Vats, D., Flegal, J. M., and, Jones, G. L. (2015). "Multivariate Output Analysis for Markov chain Monte Carlo". arXiv preprint arXiv:1512.07713.
#' 
#' 
#' @examples
#' 
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
#' @export
#'

stats.aMTM <- function(X){
   #mean squared (euclidian) jump distance
   dif <- diff(X)
   msejd <- mean(sqrt(apply(dif^2,1,sum)))
   Sigma <- var(X)
   Sigmainv <- solve(Sigma)
   msjd <-  mean(sqrt( apply(dif, 1, function(row) row %*% Sigmainv %*% row) ))
   #autocorrelation time of the mean
   SigmaP <- mcmcse::mcse.multi(X)$cov
   S <- t(chol(Sigma))
   Sinv <- solve(S)
   ACT <- Sinv %*% SigmaP %*% t(Sinv)
   act <- sqrt(sum(ACT^2))#frobenius of ACT
   #multivariate ESS
   ess <- mcmcse::multiESS(X)
   #output
   round(c(msejd = msejd, msjd = msjd,act = act, ess=ess),3)
}


