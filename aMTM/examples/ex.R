library(aMTM)
# Banana log-density with parameter B and a
p <- function(x, p) apply(x,1,function(x) -x[1]^2/(2*p$a^2) - 1/2*(x[2]+p$B*x[1]^2-p$B*p$a^2)^2)
# setup
set.seed(1)
N<-1e5;K<-3
B<-0.04;a<-8
# aMTM sampling with ASWAM update
mcmc <- aMTM(target=p, N=N, K=K, x0=c(0,0), parms=list(a=a,B=B), burnin=0.1)



#######################
# for aMTM


# acceptance rate (target is 0.5)
mcmc$acc.rate

# plot 1D samples
plot(mcmc$X)
# plot 2D sample
plot(as.matrix(mcmc$X))
# add the final covariances
for(k in seq(K)) mixtools::ellipse(mcmc$mu[,k], mcmc$Sig[,,k]*mcmc$lam[k],alpha = 0.1, col=k)

# estimate of the mean (true value is (0,0))
colMeans(mcmc$X)
# estimate of the variance and true value
round(cov(mcmc$X),4)
matrix(c(a^2,0,0,1+B^2*a^4*2),2,2,T)

# estimation of the mean with MC standard error
mcmcse::mcse.mat(mcmc$X)
# asymptotic covariance of the estimator
mcmcse::mcse.initseq(mcmc$X)$cov
# multivariate effective sample size
mcmcse::multiESS(mcmc$X)





#######################
# for plot.aMTM

#plot the pairs showing color
plot.aMTM(mcmc, color=T)
#plot the pairs showing jumps and color
plot.aMTM(mcmc, type='l', color=T)
#plot the marginals with colors
plot.aMTM(mcmc, pairs=F, color=T)


#######################
# for stats.aMTM
stats.aMTM(mcmc)
