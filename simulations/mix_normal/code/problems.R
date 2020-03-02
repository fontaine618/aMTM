#############################################
# Definition of the target density
set.seed(1)
# dimension
d <- 5
# number of mixture components
M <- 3
# weight of the components
w.t <- rep(1/M,M)
w.t <- w.t/sum(w.t)
# matrix of means
mu.t <- matrix(rnorm(M*d,0,0.1), M, d, T)
# array of covariances
Sigma.t <- array(0, dim = c(d,d,M))
for(m in seq(M)){
   S <- matrix(rnorm(d*d), d, d, T)/10
   Sigma.t[,,m] <- S %*% t(S)
}
# precomputation to speed up the evaluations
S <- sapply(seq(M), function(i)chol(Sigma.t[,,i]), simplify = 'array')
Sinv <- sapply(seq(M), function(i)solve(Sigma.t[,,i]), simplify = 'array')
denum <- sapply(seq(M), function(i)sqrt(det(Sigma.t[,,i]* 2*3.14159)) )
B <- 10^max(round(log(denum, 10)))
denum <- denum/B
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

# check that it actually works
p(t(parms$mu),parms)
# precompile for faster evaluation
target <- compiler::cmpfun(p)


# prepare data and function
data = list(
   parms = parms,
   target = p
)
fun = function(data, job, ...){
   list()
}


K = 3
d = dim(parms$mu)[1]
x0 = rep(0,d)
data$target(matrix(x0,K,d,T),data$parms)
data$target(t(data$parms$mu),data$parms)
