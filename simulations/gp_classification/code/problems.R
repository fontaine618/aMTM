#############################################
# Generating the problem 

p = 3
l2 = rnorm(p)^2
theta = log(l2)

n = 10
X = matrix(runif(n*p),n,p)

K = matrix(NA, n, n)

for(i in seq(n)){
   for(j in seq(n)){
      d = X[i,] - X[j,]
      K[i,j] = exp(- 0.5 * sum(d*d/l2))
   }
}

Xd = array(NA, dim=c(n, n, p))
for(i in seq(n)){
   for(j in seq(n)){
      Xd[i,j,] = (X[i,] - X[j,])^2
   }
}


 
f = mvtnorm::rmvnorm(1, rep(0,n), sigma=K)
w1 = exp(-f)
wm1= exp(f)
prob1 = w1 / (w1 + wm1)
y = rbinom(n, 1, prob1)



#############################################
# Define the target


K = 2
nimp = 7
theta = matrix(rnorm(p*K)^2, K, p)

parms = list(n=n, p=p, y=y, Xd=Xd, nimp=nimp)

logp = function(theta, parms){
   K = nrow(theta)
   l2 = exp(theta)
   llk = sapply(seq(K), function(k){
      K = exp(-0.5*apply(parms$Xd, 1:2, function(x) sum(x/l2[k,])))
      f = mvtnorm::rmvnorm(parms$nimp, rep(0,n), sigma=K)
      w1 = exp(-f)
      wm1= exp(f)
      py1 = w1 / (w1 + wm1)
      pf = mvtnorm::dmvnorm(f, rep(0,n), sigma=K)
      qf = pf
      py = apply(parms$y * py1 + (1-parms$y) * (1-py1), 1, prod)
      pyhat = mean(py * pf / qf)
   })
   prior = 0.0
}













#############################################
# Tests



library(aMTM)

mcmc = aMTM::aMTM(
   target = data$target,
   x0 = Xu + rnorm(8, 0, 0.1),
   parms = data$parms,
   K=3,
   N=1e4,
   global=F,
   local=F,
   scale=T,
   accrate=0.3,
   proposal=2,
   gamma=0.9,
   adapt=2
)

plot.aMTM(mcmc, vars = c(1,2), type='b')


points(as.matrix(mcmc$X[,1:2]), col=5)
points(as.matrix(mcmc$X[,3:4]), col=2)
points(as.matrix(mcmc$X[,5:6]), col=3)
points(as.matrix(mcmc$X[,7:8]), col=4)


mcmc$sel.prop
mcmc$acc.rate
