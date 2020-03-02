library(batchtools)
library(snow)
setwd("~/git/aMTM/simulations/sensor_network")

reg = makeExperimentRegistry(file.dir = NA, seed = 1)
reg$cluster.functions = makeClusterFunctionsSocket(ncpus = 15)
# add problem
source('code/problems.R')
addProblem(name = "mix_normal", data = data, fun = fun)
# add algorithms
source("code/algorithms.R")
addAlgorithm(name = "aMTM", fun = aMTM.wrapper)
# add experiments
source("code/experiments.R")
addExperiments(pdes, ades, repls = 5)

summarizeExperiments(by = c("problem", "algorithm"))

testJob(1, external = FALSE)
testJob(1, external = TRUE)

submitJobs()
getStatus()


findErrors()
getErrorMessages()

loadResult(1)

id = function(res) res$stats
results = unwrap(reduceResultsDataTable(fun = id))

pars = unwrap(getJobPars())
tab = ijoin(pars, results)
class(tab)
head(tab)


tab[,list(mean.act = mean(act)),by=(K)]
