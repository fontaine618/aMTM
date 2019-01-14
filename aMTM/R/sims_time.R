####################################################
### EXPERIMENTS COMPARING COMPLEXITY MEASURES
###
### We compare CPU time to Number of target evaluations
### for Metropolis, ASWAM (i.e. K=1), MTM and aMTM (K=2,3,5,10,20)
### given a cheap target and an expansive target.
### We do not really care about the actual output or efficiency,
### so the tuning is irrelevent here.
###
### Last update : 13-01-2019
####################################################

### Global parameters
# dimension
d <- 5
# MC smaple size
N <- 1e2
# number of modes
M <- 100

### Mixture density
set.seed(1)
w.t <- rep(1/M,M)
w.t <- w.t/sum(w.t)
id <- sample.int(M, N, replace=T, prob = w.t )
mu.t = array(rnorm(d*M, 0,3), dim = c(M,d))
Sigma.t = array(0, dim = c(d,d,M))
for(m in seq(M)){
   S <- matrix(rnorm(d^2,0,2),d)
   Sigma.t[,,m] <- S%*%t(S)+0.01*diag(d)
}
S <- sapply(seq(M), function(i)chol(Sigma.t[,,i]), simplify = 'array')
Sinv <- sapply(seq(M), function(i)solve(Sigma.t[,,i]), simplify = 'array')
denum <- sapply(seq(M), function(i)sqrt(det(Sigma.t[,,i])) * 2*3.14159)
parms <- list(M=M,mu=t(mu.t),Sinv=Sinv,w=w.t,denum=denum)
p <- function(x,pa){
   x <- as.matrix(x)
   apply(x,1,function(x){
      log(sum(sapply(seq(pa$M), function(i){
         pa$w[i] * exp(-0.5*(x-pa$mu[,i]) %*% pa$Sinv[,,i] %*% (x-pa$mu[,i]))/pa$denum[i]
      })))})
}

parms100 <- parms
p100 <- p
### average CPU time of one evaluation
NN <- 100000
system.time(p(matrix(0,NN,d), parms))/NN
#M=100
#      user    system   elapsed 
# 0.0028184 0.0000023 0.0028492
#M=1
#      user    system   elapsed 
# 0.0001142 0.0000007 0.0001175
# about 25x longer to compute.


### IID sample
x <- matrix(rnorm(N*d,0,1),N,d)
x <- sapply(seq(nrow(x)),function(i) {
   mu.t[id[i],] + t(S[,,id[i]]) %*% x[i,]
})
#pairs(t(x))
ranges <- apply(t(x),2,range)



### Timings
library(aMTM)
B <- 10
Ks <- seq(10)
times <- array(NA, dim = c(length(Ks),2,2, B), 
               dimnames = list(Ks, c('NoAdapt', 'Adapt'), c(1,100), seq(B)))
# first index is numbers of tries, 2nd is adapt/no adapt, 3rd is replication
for(b in seq(B)){
   for(k in seq(length(Ks))){
      times[k,2,1,b] <- aMTM(target=p1, N=N, K=Ks[k], x0=rep(0,5), parms=parms1)$time
      times[k,1,1,b] <- aMTM(target=p1, N=N, K=Ks[k], x0=rep(0,5), parms=parms1, adapt=0)$time
      times[k,2,2,b] <- aMTM(target=p100, N=N, K=Ks[k], x0=rep(0,5), parms=parms100)$time
      times[k,1,2,b] <- aMTM(target=p100, N=N, K=Ks[k], x0=rep(0,5), parms=parms100, adapt=0)$time
   }
}
times.mean <- apply(times, 1:3, mean)
times.min <- apply(times, 1:3, min)
times.max <- apply(times, 1:3, max)
matplot(2*Ks, times.mean[,,1]/0.0001175/N, type = 'b')
abline(0,1)
matplot(2*Ks, times.mean[,,2]/0.0028148/N, type = 'b')
abline(0,1)
