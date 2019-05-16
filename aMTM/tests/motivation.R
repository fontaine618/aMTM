library(aMTM)
library(RColorBrewer)
library(SamplerCompare)
library(car)
setwd("C:/Users/Simon/Dropbox/simon_memoire/")
###############################################################################
## Bimodal example (d=2)
###############################################################################
N <- 10000
set.seed(10)
pdf(file='tex/figures/aMTM/motivation2.pdf', width =6, height=8)
par(mfrow=c(3,2),mar=c(4,4,2,2), cex = 0.6)
#############################
## Target density
mu1 <- c(20,0)
mu2 <- c(0,8)
Sig1 <- matrix(c(9,0,0,1),2,2,T)
Sig2 <- matrix(c(1,0,0,9),2,2,T)
S1 <- solve(Sig1)
S2 <- solve(Sig2)
d1 <- (det(2*pi*Sig1))^-0.5
d2 <- (det(2*pi*Sig2))^-0.5
w1 <- 0.3
w2 <- 1-w1
pi <- function(x){
   w1 * d1 * exp(- (x-mu1)%*%S1%*%(x-mu1) /2)+
   w2 * d2 * exp(- (x-mu2)%*%S2%*%(x-mu2) /2)
}
logpi <- function(x){
   log(
      w1 * d1 * exp(- t(x-mu1)%*%S1%*%(x-mu1) /2)+
         w2 * d2 * exp(- t(x-mu2)%*%S2%*%(x-mu2) /2)
   )
}
grid <- expand.grid(x=seq(-5,30,0.1), y= seq(-5,20,0.1))
Z <- apply(grid, 1, pi)
z <- matrix(Z,351,251)
#contour(seq(-5,30,0.1),seq(-5,15,0.1),z,nlevels=5)
#abline(h=0,lty=3)
#abline(v=0,lty=3)
image(x=seq(-5,30,0.1), y=seq(-5,20,0.1), z=z,
      col = grey(seq(1,0,-0.1)), bty='o',
      xlab = expression(x[1]),ylab = expression(x[2]),
      main = '(a) Densité cible')
axis(2,labels=NA,lwd.ticks=0)
axis(3,labels=NA,lwd.ticks=0)
abline(h=0,lty=3)
abline(v=0,lty=3)
#############################
## IID sample
K <- sample.int(2, N, T, c(w1,w2))
X <- t(sapply(K, function(k){
   z <- rnorm(2)
   if(k==1) z[1]<-z[1]*3 + 20
   if(k==2) z[2]<-z[2]*3 + 8
   z
}))
X.iid <- X
plot(X, xlim=c(-5,30),ylim=c(-5,20),
     xlab = expression(x[1]),ylab = expression(x[2]),
     main = '(b) Échantillon i.i.d.' ,xaxs='i',yaxs='i')
abline(h=0,lty=3)
abline(v=0,lty=3)
prop <- round(sum(X[,1]>5)/N,3)
legend("topright", legend =
          c(as.expression(substitute(Pr(x[1]>5) == pr,
                                     list(pr = as.name(prop))))), bty="n") 
#############################
## AM sampler
d=2

target <- make.dist(2, '2D Bimodal Gaussian Mixtue', '2D Bimodal Gaussian Mixture',
                    logpi)

out <- adaptive.metropolis.sample(target, c(0,10),N,tuning=1)
X.AM1 <- out$X
plot(out$X, xlim=c(-5,30),ylim=c(-5,20),
     xlab = expression(x[1]),ylab = expression(x[2]),
     main = '(c) Échantillon AM' ,xaxs='i',yaxs='i')
abline(h=0,lty=3)
abline(v=0,lty=3)
ellipse(apply(out$X,2,mean), out$sample.cov,radius=1.96*(2.38)^2/d)
alp <- round(1-out$reject.rate,3)
prop <- round(sum(out$X[,1]>5)/N,3)
legend("topright", legend =
          c(as.expression(substitute(alpha == a,
                                     list(a = as.name(alp)))),
            as.expression(substitute(Pr(x[1]>5) == pr,
                                     list(pr = as.name(prop))))), bty="n") 

out <- adaptive.metropolis.sample(target, c(20,0),N,tuning=20)
X.AM2 <- out$X
plot(out$X, xlim=c(-5,30),ylim=c(-5,20),
     xlab = expression(x[1]),ylab = expression(x[2]),
     main = '(d) Échantillon AM' ,xaxs='i',yaxs='i')
abline(h=0,lty=3)
abline(v=0,lty=3)
ellipse(apply(out$X,2,mean), out$sample.cov,radius=1.96*(2.38)^2/d)
alp <- round(1-out$reject.rate,3)
prop <- round(sum(out$X[,1]>5)/N,3)
legend("topright", legend =
c(as.expression(substitute(alpha == a,
                           list(a = as.name(alp)))),
  as.expression(substitute(Pr(x[1]>5) == pr,
                           list(pr = as.name(prop))))), bty="n") 
#############################
## MTM sampler
source("R/tests/mtm_sampler.R")
logtarget=logpi

x0<-c(0,10)
K<-3
d<-2
N<-N
cov<-array(c(40,0,0,4,
             100,0,0,100,
             4,0,0,40), dim=c(d,d,K))
out <- mtm(logtarget, x0, K, cov, N, d)
X.MTM3 <- out$X
plot(out$X, xlim=c(-5,30),ylim=c(-5,20),col=out$Sel,
     xlab = expression(x[1]),ylab = expression(x[2]),
     main = '(e) Échantillon MTM' ,xaxs='i',yaxs='i')
abline(h=0,lty=3)
abline(v=0,lty=3)
for(k in seq(K)){
   ellipse(apply(out$X,2,mean), cov[,,k],radius=1.96,col=k)
}
alp <- round(mean(out$Acc),3)
prop <- round(sum(out$X[,1]>5)/N,3)
legend("topright", legend =
          c(as.expression(substitute(alpha == a,
                                     list(a = as.name(alp)))),
            as.expression(substitute(Pr(x[1]>5) == pr,
                                     list(pr = as.name(prop))))), bty="n") 


x0<-c(0,10)
K<-2
d<-2
N<-N
cov<-array(c(140,-50,-50,30,
             8,0,0,8), dim=c(d,d,K))
out <- mtm(logtarget, x0, K, cov, N, d)

X.MTM2 <- out$X
plot(out$X, xlim=c(-5,30),ylim=c(-5,20),col=out$Sel,
     xlab = expression(x[1]),ylab = expression(x[2]),
     main = '(f) Échantillon MTM' ,xaxs='i',yaxs='i')
abline(h=0,lty=3)
abline(v=0,lty=3)
for(k in seq(K)){
   ellipse(apply(out$X,2,mean), cov[,,k],radius=1.96,col=k)
}
alp <- round(mean(out$Acc),3)
prop <- round(sum(out$X[,1]>5)/N,3)
legend("topright", legend =
          c(as.expression(substitute(alpha == a,
                                     list(a = as.name(alp)))),
            as.expression(substitute(Pr(x[1]>5) == pr,
                                     list(pr = as.name(prop))))), bty="n") 
dev.off()

#############################
## aMTM package application
mu1 <- c(20,0)
mu2 <- c(0,8)
Sig1 <- matrix(c(9,0,0,1),2,2,T)
Sig2 <- matrix(c(1,0,0,9),2,2,T)
S1 <- solve(Sig1)
S2 <- solve(Sig2)
d1 <- (4*3.1416*det(Sig1))^-0.5
d2 <- (4*3.1416*det(Sig2))^-0.5
w1 <- 0.3
w2 <- 1-w1
N <- 10000


logp <- function(X, p){
   if(is.vector(X)) X <-t(X)
   X <- as.matrix(X)
   apply(X, 1, function(x){
     log(p$w1 * p$d1 * exp(- t(x-p$mu1)%*%p$S1%*%(x-p$mu1) /2)+
         p$w2 * p$d2 * exp(- t(x-p$mu2)%*%p$S2%*%(x-p$mu2) /2)
     )})
   }
#logp <- compiler::cmpfun(logp)
p <- list(w1=w1,w2=w2,d1=d1,d2=d2,mu1=mu1,mu2=mu2,S1=S1,S2=S2)
#X <- matrix(c(0,0,20,5,10,0),3,2)
K<-2; d<-2
sig0 = array(0, dim = c(d,d,K))
for(k in seq(K)) sig0[,,k] <- diag(d)*10^k


library(aMTM)

set.seed(1)
mcmc <- aMTM(target=logp, N=1e4, K=K, x0=c(0,10), parms=p,
             sig0 = sig0, adapt = 'ASWAM', local = T, accrate = 0.5)



mean(mcmc$X[,1]>5)
mcmc$sel.prop
mcmc$acc.rate
mcmc$time



setwd("C:/Users/Simon/Dropbox/simon_memoire/tex/figures/sims")
pdf(file='motiv_trace_plots.pdf', width =6, height=4)
par(cex=0.6)
plot.aMTM(mcmc, pairs=F)
dev.off()


setwd("C:/Users/Simon/Dropbox/simon_memoire/tex/figures/sims")
pdf(file='motiv_pairs_plots.pdf', width =6, height=6)
par(cex=0.6)
plot.aMTM(mcmc, type='b')
dev.off()

stats.aMTM(mcmc$X)
coda::autocorr.plot(mcmc$X)
coda::geweke.plot(mcmc$X)

(diag <- coda::geweke.diag(mcmc$X))
# Fraction in 1st window = 0.1
# Fraction in 2nd window = 0.5 
# 
#    var1    var2 
# -0.1348 -0.2013
pnorm(abs(diag$z),lower.tail=FALSE)*2
#      var1      var2 
# 0.8927779 0.8404725 

coda::rejectionRate(mcmc$X)
coda::HPDinterval(mcmc$X)
#          lower    upper
# var1 -2.279166 23.46645
# var2 -1.637213 12.95471
coda::effectiveSize(mcmc$X)
#     var1      var2 
# 79.17247 113.43820 

stats.aMTM(mcmc$X)
# msejd    msjd     act     ess 
# 1.109   0.325  66.332 409.908

setwd("C:/Users/Simon/Dropbox/simon_memoire/tex/figures/sims")
pdf(file='motiv_acfplot.pdf', width =6, height=3)
par(mfrow=c(1,1), cex = 0.6)
coda::acfplot(mcmc$X, lag.max=200, aspect = '')
dev.off()

set.seed(2)
mcmc2 <- aMTM(target=logp, N=1e4, K=K, x0=c(0,10), parms=p,
             sig0 = sig0, adapt = 'ASWAM', local = T, accrate = 0.5)
mcmc.list <- coda::mcmc.list(mcmc$X, mcmc2$X)
coda::gelman.plot(mcmc.list)
coda::gelman.diag(mcmc.list)

mcmclist <- lapply(seq(10), function(i){
   set.seed(i)
   mcmclist[[i]] <- aMTM(target=logp,N=1e4,K=K,x0=c(0,10), 
                         parms=p,sig0=sig0,adapt='ASWAM',
                         local=T,accrate=0.5)$X
})
coda::gelman.diag(mcmclist)
# Potential scale reduction factors:
#    
#      Point est. Upper C.I.
# [1,]       1.02       1.03
# [2,]       1.01       1.02
# 
# Multivariate psrf
# 
# 1.02

mcmcse::mcse.multi(mcmc$X)$est
#     var1     var2 
# 5.834973 5.781332 
mcmcse::mcse.multi(mcmc$X)$cov/N
#            [,1]       [,2]
# [1,]  0.5686868 -0.2347880
# [2,] -0.2347880  0.1035699
mcmcse::mcse.mat(mcmc$X)
