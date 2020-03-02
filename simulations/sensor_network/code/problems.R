#############################################
# Generating the problem 2, 22, 6, 25, 26 38!

n_sensors = 6
n_known = 2
n_unknown = 4
known_sensors = c(F,F,F,F,T,T)
set.seed(6)
X = round(matrix(
   c(
      0.57,0.91,
      0.10,0.37,
      0.24,0.14,
      0.85,0.04,
      0.50,0.30,
      0.30,0.70
      ),
   nrow=n_sensors,
   ncol=2, T
),2)
y = as.matrix(dist(X, diag=F, upper=F))

probs = exp(-y^2 / (2*0.5^2))
w = apply(probs, 1:2, function(p) rbinom(1,1,p))
w[upper.tri(w)] = 0

w = diag(n_sensors)
w[c(4,5,6), 1] = 1
w[c(3,4), 2] = 1
w[c(5,6), 3] = 1


y = y*w
y[y==0]=NA
y[upper.tri(y)] = NA

segs = data.frame(i=NA, j=NA, x0=NA, y0=NA, x1=NA, y1=NA)


for(i in seq(2,n_sensors)){
   for(j in seq(1, i-1)){
      if(w[i,j] == 1){
         segs = rbind(segs, c(i,j,X[i,1],X[i,2],X[j,1],X[j,2]))
      } 
   }
}


# Define the target
Xk = as.vector(t(X[known_sensors,]))
Xu = as.vector(t(X[!known_sensors,]))

parms=list(Xk=Xk, y=y, w=w, sigy=0.02, sigw=0.3, sigx=1000, 
           ns=n_sensors, nu=n_unknown, nk=n_known)

logp = function(x, parms){
   K = length(x)/(2*parms$nu)
   xu= x
   x = cbind(x, matrix(parms$Xk, K, parms$nk*2, byrow=T))
   X = array(x, dim=c(K, 2, parms$ns))
   
   out = sapply(seq(K), function(k){
      d = as.matrix(dist(t(X[k,,]), diag=T, upper=T))
      llky = -sum((parms$y - d)^2 / (2*parms$sigy^2), na.rm=T)
      llkw = -sum(parms$w*d^2 / (4*parms$sigw^2), na.rm=T)
      llkw1 = sum((1-parms$w) * log(1-exp(-d^2 / (2*parms$sigw^2))), na.rm=T)
      prior = -exp(sum((xu[k, ]-0.5)^2) / (2*parms$sigx))
      #prior = -sum(abs(xu[k, ]) > 1) * 1e+10
      return(llky + llkw + llkw1 + prior)
   })
   return(out)
}


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

K = 3
x = Xu + rnorm(2*n_unknown*K, 0, 0.5)
x = matrix(x, K, 2*n_unknown, byrow=T)
data$target(x,data$parms)







plot_target = function(){
   
   par(mfrow=c(1,1))
   
   plot(X, xlim=c(-1,1.5), ylim=c(-0.5,1.5), 
        col=known_sensors+1, pch=19)
   text(x=X[,1]+0.05, y=X[,2]+0.05, labels=1:n_sensors)
   segments(segs$x0, segs$y0, segs$x1, segs$y1)
   abline(h=0, lty=3)
   abline(v=0, lty=3)
}



plot_target()






library(aMTM)

K = 5

log_seq = function(start, end, n){
   l = seq(log(start), log(end), length.out = n)
   exp(l)
}

sig0 = sapply(log_seq(1, 0.01, K), function(sig2){
   diag(rep(1,n_unknown*2)^2) * sig2
}, simplify="array")



N = 2e5
x0 = Xu + rnorm(2*n_unknown, 0, 0.02)
x0 = rnorm(2*n_unknown, 0.5, 0.2)

mcmc = aMTM::aMTM(
   target=data$target,
   x0=x0,
   parms=data$parms,
   sig0=sig0,
   K=K,
   N=N,
   global=T,
   local=F,
   scale=F,
   accrate=0.3,
   proposal=3,
   gamma=0.7,
   adapt=2,
   burnin=0,
   weight=-1
)

ids = seq(1,N,by=100)
plot_target()
points(as.matrix(mcmc$X[ids,1:2]), col=1)
points(as.matrix(mcmc$X[ids,3:4]), col=2)
points(as.matrix(mcmc$X[ids,5:6]), col=3)
points(as.matrix(mcmc$X[ids,7:8]), col=4)


mcmc$sel.prop
mcmc$acc.rate

round(mcmc$Sig, 5)


data$target(mcmc$X[1:10, ],data$parms)


plot.aMTM(mcmc, vars=c(1,2), type="b")
plot.aMTM(mcmc, vars=c(3,4), type="b")
plot.aMTM(mcmc, vars=c(5,6), type="b")
plot.aMTM(mcmc, vars=c(7,8), type="b")


