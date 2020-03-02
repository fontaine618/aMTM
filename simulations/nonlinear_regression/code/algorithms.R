aMTM.wrapper = function(data, job, instance, ...){
   mcmc = aMTM::aMTM(
      target = data$target,
      x0 = rep(0,dim(data$parms$mu)[1]),
      parms = data$parms, ...
   )
   list(stats = aMTM::stats.aMTM(mcmc$X))
}
