burnin <- 1000*40
N <- 19000*40
thin <- 1
nchains <- 4
readChains <- function(filename) {
  XX <- array(scan(filename), dim=c(4, N, nchains))
  XX <- aperm(XX, c(2,1,3))
  X <- list()
  library(coda)
  transform <- function(X) cbind(plogis(X[,1:3]), X[,4])
  for(k in seq_len(dim(XX)[3]))
    X[[k]] <- mcmc(transform(XX[,,k]), start=(burnin*thin), thin=thin)
  xx <- do.call(mcmc.list, X)
  return(xx)
}
fnames <- list(RAPTOR="chains-RAPTOR.dat", AM="chains-AM.dat")
LOH <- lapply(fnames, readChains)

LOH$AM <- LOH$AM[c(1,3)] #remove bad chains
lapply(LOH, function(xx) {
  ans <- t(apply(as.matrix(xx), 2, function(x) c(mean=mean(x), sd=sd(x))))
  rownames(ans) <- c("eta", "pi1", "pi2", "gamma")
  return(ans)
})

mean.abs.acf <- array(NA, dim=c(4, 2, 4))
for(row in 1:4) for(col in 1:2) for(ch in seq_along(LOH[[col]])) {
  z <- LOH[[col]][[ch]][,row]
  mean.abs.acf[row, col, ch] <- mean(abs(acf(z, lag.max=40*40, plot=FALSE)$acf))
}
apply(mean.abs.acf, 1:2, mean, na.rm=TRUE)

LOH <- as.matrix(LOH$RAPTOR)
set.seed(1234)
LOH <- LOH[sample(20000),]
save(LOH, file="LOH.RData")
