#############################################
# Generating the problem 
setwd("~/git/aMTM/simulations/loh/code")
df = read.table("../data/BarrettsLOH.dat", header=F)
colnames(df) = c("x", "n")


parms = list(
   x=df$x, 
   n=df$n,
   gam_bound=30,
   dbetabinom=rmutil::dbetabinom
)



logp = function(theta, parms){
   K = nrow(theta)
   eta = pmin(pmax(theta[,1], 1.0e-12),1-1.0e-12)
   pi1 = pmin(pmax(theta[,2], 1.0e-12),1-1.0e-12)
   pi2 = pmin(pmax(theta[,3], 1.0e-12),1-1.0e-12)
   gam = theta[,4]
   w = exp(gam) /(2*(1+exp(gam))) 
   cmp1 = sapply(pi1, function(pi) 
      dbinom(parms$x, parms$n, pi)
      )
   cmp2 = mapply(function(pi, s) 
      parms$dbetabinom(parms$x, parms$n, pi, 1/s),
      pi2, w
      )
   eta_mat = matrix(eta, length(parms$x), K, T)
   f = eta_mat * cmp1 + (1-eta_mat) * cmp2
   llk = apply(log(f), 2, sum)
   prior = rep(0, K)
   prior = prior - (eta < 1e-10)
   prior = prior - (eta > 1 - 1e-10)
   prior = prior - (pi1 < 1e-10)
   prior = prior - (pi1 > 1 - 1e-10)
   prior = prior - (pi2 < 1e-10)
   prior = prior - (pi2 > 1 - 1e-10)
   prior = prior - (abs(gam) > parms$gam_bound)
   return(llk + prior*1e5 + 100)
}


p = 4
K = 7
theta = matrix(
   c(runif(K), runif(K), runif(7), runif(K, -30, 30)),
   K, p, F 
)

logp(theta, parms)



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


plotContours = function(){
   par(mfrow=c(6,6), oma=rep(0,4), mar=rep(0,4))
   
   for(eta in c(0.01,0.1,0.25,0.5,0.8,0.99)){
      for(gam in c(-25,-15,-5,5,15,25)){
         theta = as.matrix(expand.grid(
            eta = eta,
            pi1 = seq(0.001, 0.999, length.out=100),
            pi2 = seq(0.001, 0.999, length.out=100),
            gamma = gam
         ))
   
         vals = matrix(logp(theta, parms), 100, 100)
   
         contour(
            x=seq(0.001, 0.999, length.out=100),
            y=seq(0.001, 0.999, length.out=100),
            z=vals,
            levels=seq(-50, 10, 10)
         )
      }
   }
}

# plotContours()




#############################################

theta0 = c(0.078, 0.83, 0.23, -18)
theta0 = c(0.897, 0.229, 0.714, 15.661)

theta0 = c(runif(3), rnorm(1,0,10))

K = 5

log_seq = function(start, end, n){
   l = seq(log(start), log(end), length.out = n)
   exp(l)
}

sig0 = sapply(log_seq(1, 0.001, K), function(sig2){
   diag(c(0.2,0.2,0.2,10)^2) * sig2
}, simplify="array")

N = 1e5

mcmc = aMTM::aMTM(
   target = data$target,
   x0 = theta0,
   parms = data$parms,
   sig0=sig0,
   K=K,
   N=N,
   global=T,
   local=F,
   scale=F,
   accrate=0.3,
   proposal=2,
   gamma=0.5,
   adapt=2,
   burnin=0,
   weight=-1
)

aMTM::plot.aMTM(mcmc, vars = 1:4, type='b')
aMTM::plot.aMTM(mcmc, vars = 1:4, pairs=F)

mcmc$sel.prop
mcmc$acc.rate

mcmc$Sig
theta0

ids = seq(1,N,by=100)
plot(as.matrix(mcmc$X[ids,2:3]), col=1, xlim=0:1, ylim=0:1)



