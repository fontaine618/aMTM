#WRONG!, maybe redu with arguyment f= function of choice, default is Id
# use frobenius norm on dimension of f
ACT <- function(x, max.lag = min(100,nrow(x)-2)){
   n <- nrow(x)
   m <- apply(x,2,mean)
   xc <- t(t(x) - m)
   gamma <- sapply(seq(max.lag), function(i){
      sum( apply( xc[1:(n-i),] * xc[(i+1):n,],2,sum) ) / (n-i)
   })
   gam0 <- sum( apply( xc[1:n,] * xc[1:n,],1,sum) ) / (n)
   1+2*sum(gamma)/gam0
}

MSEJD <- function(x){
   mean(sqrt(apply(diff(x)^2,1,sum)))
}

plot.pairs <- function(mcmc, mu.t, Sigma.t, burnin, ranges, dim, tmp){
   d <- ncol(mu.t)
   N <- nrow(mcmc$X)
   M <- nrow(mu.t)
   K <- mcmc$K
   if(missing(dim)) dim <- seq(d)
   dim <- intersect(dim,seq(d))
   if(length(dim)>6)dim<-dim[1:6]
   par(mfrow = c(length(dim),length(dim)), mar= rep(0.1,4),
       oma = c(3,3,4,3))
   if(missing(ranges)) ranges <- apply(mcmc$X,2,range)
   plots <- expand.grid(j=dim,i=dim)



   if(burnin<=1){
      ids <- seq(N*burnin+1,N)
   }else{
      ids <- seq(burnin+1,N)
   }
   for(k in seq(nrow(plots))){
      i <- plots[k,2]
      j <- plots[k,1]
      xlim <- ranges[,j]
      ylim <- ranges[,i]
      plot(NA, xlim = xlim, ylim=ylim,
                    yaxt = ifelse(j==1,'s','n'), xaxt = ifelse(i==max(dim),'s','n'))
      if(i==j){
         text(x=ranges[1,i],y=ranges[2,i],pos=1,labels=i,cex=2)
         if(i==1){
            dat <- paste(c('time','msejd','act','dist.mean','miss.mode','dist.TV*M','acc.rate'),
                         round(tmp[7:13],3), sep=' = ')
            mtext(text = paste(dat, collapse = '     '),
                  side =3, outer=T, cex = 0.8)
         }
      }else{
         points(x=mcmc$X[ids,j],y=mcmc$X[ids,i],pch=1,col=mcmc$sel[ids])
         for(k in seq(K)){
            mu <- apply(mcmc$X[ids,c(j,i)],2,median)
            #sig <-  mcmc$Sig[N,k,c(j,i),c(j,i)]
            sig <-  mcmc$Sig[k,c(j,i),c(j,i)]
            ellipse(mu=mu, sigma=sig, alpha = .05,
                    npoints = 250, col=k, lty = 3, lwd =2)
         }
         for(m in seq(M)){
            mu <- mu.t[m,c(j,i)]
            sig <- Sigma.t[c(j,i),c(j,i),m]
            ellipse(mu=mu, sigma=sig, alpha = .01,
                    npoints = 250, col=1, lty = 1,lwd=2)
         }
      }
   }
   dat <- paste(c('alg','common','gamma','target.acc','K'),
                round(tmp[2:6],3), sep=' = ')
   title(main=paste(dat, collapse = '     '), outer=T)
}

plot.rep <- function(x, y, ...){
   y <- y[!is.na(x)]
   x <- x[!is.na(x)]
   x.vals <- unique(x)
   n <- length(x.vals)
   m <- rep(NA,n)
   se <- rep(NA,n)
   for(i in seq(n)){
      m[i] <- mean(y[x==x.vals[i]])
      se[i] <- sd(y[x==x.vals[i]])/sqrt(length(y[x==x.vals[i]]))
   }
   plot(NA, xlim = range(x.vals, na.rm=T), ylim = range(y, na.rm=T), ...)
   lines(x.vals,m,lty=1)
   lines(x.vals,m+se,lty=3)
   lines(x.vals,m-se,lty=3)
}

mix.compare <- function(mcmc,ids,parms,x,mlda){
   X <- mcmc$X[ids,]
   M <- parms$M
   #mean sqaured expected jump distance
   msejd <- MSEJD(X)
   #autocorrelation time of the mean
   act <- ACT(X)
   #bias of the mean in mahalanobis distance
   m.exp <- apply(X,2,mean)
   m.true <- apply(parms$mu,1,mean)
   Sigma.total <- var(t(x))
   mean.dist <- mahalanobis(m.exp, m.true, Sigma.total)
   #mlda
   pred <- predict(mlda, X)
   prop <- table(pred$class)/nrow(X)
   miss.mode <- sum(abs(prop - mlda$counts/mlda$N))/2
   #TV distance
   dist.TV <- max(abs(prop - mlda$counts/mlda$N))
   #acceptation rate
   acc <- mean( mcmc$acc[ids])
   #output
   c(time=mcmc$time,
     msejd = msejd,
     act = act,
     dist.mean = mean.dist,
     miss.mode = miss.mode * M,
     dist.TV = dist.TV,
     acc.rate = acc)
}
