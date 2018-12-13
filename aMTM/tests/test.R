d <- 3
M <- 3

N <- 10000
burnin <- 0.5

# MIXTURE OF NORMALS
set.seed(4)
x <- matrix(rnorm(N*d,0,1),N,d)
w.t <- rep(1/M,M)
w.t <- w.t/sum(w.t)
id <- sample.int(M, N, replace=T, prob = w.t )
mu.t = array(rnorm(d*M, 0,10), dim = c(M,d))
Sigma.t = array(0, dim = c(d,d,M))
for(m in seq(M)){
   S <- matrix(rnorm(d^2,0,2),d)
   Sigma.t[,,m] <- S%*%t(S)+1*diag(d)
}
mu.t <- round(mu.t,0)
Sigma.t <- round(Sigma.t,0)
S <- sapply(seq(M), function(i)chol(Sigma.t[,,i]), simplify = 'array')
Sinv <- sapply(seq(M), function(i)solve(Sigma.t[,,i]), simplify = 'array')
denum <- sapply(seq(M), function(i)sqrt(det(Sigma.t[,,i])) * 2*3.14159)
x <- sapply(seq(nrow(x)),function(i) {
   mu.t[id[i],] + t(S[,,id[i]]) %*% x[i,]
})
ranges <- apply(t(x),2,range)
plot(t(x), xlim = ranges[,1], ylim = ranges[,2])

#target function
parms <- list(M=M,mu=t(mu.t),Sinv=Sinv,w=w.t,denum=denum)
p <- function(x,pa){
   x<- as.vector(x)
   log(sum(sapply(seq(pa$M), function(i){
      pa$w[i] * exp(-0.5*(x-pa$mu[,i]) %*% pa$Sinv[,,i] %*% (x-pa$mu[,i]))/pa$denum[i]
   })))
}
p <- compiler::cmpfun(p)
#initialize
#par(mfrow=c(2,2))
K <- 3
x0 = rnorm(d, 0,5)
sig0 = array(0, dim = c(d,d,K))
for(k in seq(K)){
   S <- matrix(rnorm(d^2,0,20)/2^k,d)
   sig0[,,k] <- S%*%t(S)+diag(d)*0.1
}
mu0 <- array(rnorm(d*K, 0,10), dim = c(d,K))
lam0 <- array(2.38^2/d, dim = c(K))
#c++ call
library(aMTM)
set.seed(1)
mcmc <- aMTM(target = p, N = N, K = K,
             x0 = x0, sig0 = sig0, mu0 = mu0, lam0 = lam0,
             adapt = 2, global = T, scale = T, local = T,
             proposal = 0, accrate = 0.35, gamma = 0.65,
             parms = parms, beta = 0)
X <- mcmc$X[seq(burnin*N,N),]

par(mfrow=c(2,2),mar=c(3,3,3,3))
vars <- c(1,2)
plot(X[,vars], xlim = ranges[,vars[1]], ylim = ranges[,vars[2]], col=mcmc$sel+1,
     xlab = vars[1], ylab= vars[2])
for(k in seq(K)){
   mixtools::ellipse(mu=apply(X[,vars],2,mean), sigma=mcmc$Sig[vars,vars,k]*mcmc$lam[k], alpha = .01,
                     npoints = 250, col=k, lty = 1,lwd=2)
}
vars <- c(3,2)
plot(X[,vars], xlim = ranges[,vars[1]], ylim = ranges[,vars[2]], col=mcmc$sel+1,
     xlab = vars[1], ylab= vars[2])
for(k in seq(K)){
   mixtools::ellipse(mu=apply(X[,vars],2,mean), sigma=mcmc$Sig[vars,vars,k]*mcmc$lam[k], alpha = .01,
                     npoints = 250, col=k, lty = 1,lwd=2)
}
vars <- c(1,3)
plot(X[,vars], xlim = ranges[,vars[1]], ylim = ranges[,vars[2]], col=mcmc$sel+1,
     xlab = vars[1], ylab= vars[2])
for(k in seq(K)){
   mixtools::ellipse(mu=apply(X[,vars],2,mean), sigma=mcmc$Sig[vars,vars,k]*mcmc$lam[k], alpha = .01,
                     npoints = 250, col=k, lty = 1,lwd=2)
}
plot(t(x), xlim = ranges[,1], ylim = ranges[,2])

#sapply(seq(K), function(k) mcmc$Sig[,,k] * mcmc$lam[k], simplify = 'array')
mean(mcmc$acc)
table(mcmc$sel[seq(burnin*N,N)])/(N*(1-burnin))
