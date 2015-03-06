
source("AskModelsFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("CoalescentFunctions.R")
source("PriorFunctions.R")
source("MarkovProcess.R")
source("GeneticDataSimulation.R")
source("SummaryStats.R")
source("PostAnalysis.R")

library(raster)
library(abc)
library(parallel)

###### Environmental data :

load("bio.RData")

###### Genetic data :

# load fake GeneticData 
load("dataSet.RData")

# genetic value of the ancestor
ancestor <- 200

# assuming we have the step values for each locus
steps <- c(1,2,3,4,5,4,4,4,4,4)
length(steps)

###### Model :

load("ParamList.RData")


##### 
# launch simulations
simSpatialCoal(nbSimul=10000, ParamList=ParamList, rasterStack=bio, nicheMeth = "arithmetic", GeneticData=dataSet, initialGenetValue=ancestor,
               stepValueOfLoci= steps, cores=detectCores())

# analyse the results, computes summary statistics
pcaRes <- pca4abc(GeneticData = genetic, ParamList = ParamList, distanceMethod = "DeltaMuDistance", path = paste0(getwd(),"/SimulResults"))

# performs cross validation
cv.rej <- cv4abc(param = pcaRes[[3]], sumstat=pcaRes[[2]], nval=500, tols=c(0.01, 0.1), method="rejection")
plot(cv.rej)

# performs abc analysis
result <- abc(target = pcaRes[[1]], 
              param = pcaRes[[3]], 
              sumstat = pcaRes[[2]], 
              tol=0.01, 
              method="neuralnet" )
hist(result)
