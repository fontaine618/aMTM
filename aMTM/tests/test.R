d <- 3
M <- 2

N <- 10000
burnin <- 0.1

# MIXTURE OF NORMALS
set.seed(4)
x <- matrix(rnorm(N*d,0,1),N,d)
w.t <- rep(1/M,M)
w.t <- w.t/sum(w.t)
id <- sample.int(M, N, replace=T, prob = w.t )
mu.t = array(rnorm(d*M, 0,5), dim = c(M,d))
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
#plot(t(x), xlim = ranges[,1], ylim = ranges[,2])
#pairs(t(x))
mlda <- MASS::lda(x=t(x), grouping = id)
Sigma <- var(t(x))
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
K <- 4
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
             proposal = 3, accrate = 0.30, gamma = 0.65,
             parms = parms, weight = 0, burnin=burnin)

round(mix.compare(mcmc,parms,Sigma,mlda),3)
#X <- matrix(mcmc$X,N,d)
pairs(mcmc$X, col = mcmc$sel+1)
#pairs(t(x))
mcmc$Sig
