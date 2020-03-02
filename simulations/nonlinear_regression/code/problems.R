#############################################
# Generating the problem 

x = c(1,2,3,4,5,7)
y = c(8.3,10.3,19.0,16.0,15.6,19.8)
n = length(x)
p = 2

K = 19

lam = 0.01

parms = list(x=x, y=y, n=n, p=p, lam=lam)
theta = matrix(runif(K*p), K, p) * 10

logp = function(theta, parms){
   K = nrow(theta)
   res = matrix(parms$y, parms$n, K, F) - 
      matrix(theta[,1], parms$n, K, T) * 
      (1-exp(-parms$x%*% t(theta[,2])))
   llk = - log(apply(res^2, 2, sum)) * (parms$n / 2 - 1)
   prior = (
      log((theta[,1]>-20) + 1e-20) +
      log((theta[,1]<50) + 1e-20) +
      log((theta[,2]>-2) + 1e-20) +
      log((theta[,2]<6) + 1e-20)
   )
   return(llk+prior)
}

logp(theta, parms)


theta = as.matrix(expand.grid(
   theta1 = seq(-20, 50, length.out=100),
   theta2 = seq(-2, 6, length.out=100)
))

vals = matrix(logp(theta, parms), 100, 100)

contour(
   x=seq(-20, 50, length.out=100),
   y=seq(-2, 6, length.out=100),
   z=vals
)



# precompile for faster evaluation
target <- compiler::cmpfun(logp)

# prepare data and function
data = list(
   parms = parms,
   target = logp
)
fun = function(data, job, ...){
   list()
}

#############################################
library(aMTM)


K = 5

log_seq = function(start, end, n){
   l = seq(log(start), log(end), length.out = n)
   exp(l)
}

sig0 = sapply(log_seq(1, 0.01, K), function(sig2){
   diag(c(20,2)^2) * sig2
}, simplify="array")

N = 1e5

mcmc = aMTM::aMTM(
   target = data$target,
   x0 = rep(0,p),
   parms = data$parms,
   sig0=sig0,
   K=K,
   N=N,
   global=F,
   local=T,
   scale=F,
   accrate=0.5,
   proposal=2,
   gamma=1,
   adapt=2,
   burnin=0.1,
   weight=0
)

plot.aMTM(mcmc, type='b', pairs=T)


