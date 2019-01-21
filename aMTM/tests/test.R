# to produce pdf manual
# pack <- "aMTM"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))

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
   #x <- as.vector(x)
   x <- as.matrix(x)
   apply(x,1,function(x){
   log(sum(sapply(seq(pa$M), function(i){
      pa$w[i] * exp(-0.5*(x-pa$mu[,i]) %*% pa$Sinv[,,i] %*% (x-pa$mu[,i]))/pa$denum[i]
   })))})
}
p(matrix(1:6,2,3),parms)
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
             adapt = 2, global =F, scale = F, local =T,
             proposal = 0, accrate = 0.50, gamma = 0.7,
             parms = parms, weight = 0, burnin=burnin)

round(mix.compare(mcmc,parms,Sigma,mlda),3)
#X <- matrix(mcmc$X,N,d)
pairs(as.matrix(mcmc$X), col = mcmc$sel)
mcmc$sel.prop
mcmc$lam
mcmc$Sig

pairs(t(x))

#########################################
# TEST plot.aMTM
plot.aMTM(mcmc)

stats.aMTM(mcmc$X)







#########################################
# Function for camparison with mixture target

mix.compare <- function(mcmc,parms,Sigma,mlda){
   X <- mcmc$X
   M <- parms$M
   #mean squared (euclidian) jump distance
   dif <- diff(X)
   msjd <- mean(sqrt(apply(dif^2,1,sum)))
   if(missing(Sigma)) Sigma <- var(X)
   Sigmainv <- solve(Sigma)
   msejd <-  mean(sqrt( apply(dif, 1, function(row) row %*% Sigmainv %*% row) ))
   #autocorrelation time of the mean
   SigmaP <- mcmcse::mcse.multi(X)$cov
   S <- t(chol(Sigma))
   Sinv <- solve(S)
   ACT <- Sinv %*% SigmaP %*% t(Sinv)
   act <- sqrt(sum(ACT^2))#frobenius of ACT
   #multivariate ESS
   ess <- mcmcse::multiESS(X)
   #bias of the mean in mahalanobis distance
   m.exp <- apply(X,2,mean)
   m.true <- apply(parms$mu,1,mean)
   mean.dist <- mahalanobis(m.exp, m.true, Sigma)
   if(!missing(mlda)){
      #mlda (suppress warnings beacause col names may be wrong)
      suppressWarnings(pred <- predict(mlda, X))
      prop <- table(pred$class)/nrow(X)
      miss.mode <- sum(abs(prop - mlda$counts/mlda$N))/2
      #TV distance
      dist.TV <- max(abs(prop - mlda$counts/mlda$N))
   }else{
      miss.mode=NA
      dist.TV=NA
   }
   
   #output
   c(time=mcmc$time,msjd = msjd,
     msejd = msejd,act = act, ess=ess,
     dist.mean = mean.dist, miss.mode = miss.mode * M,
     dist.TV = dist.TV,acc.rate = mcmc$acc.rate)
}
