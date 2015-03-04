
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

## Create false environmental Data
# # Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900))
# Data2 <- data.frame(BIO1=c(300,400),BIO12=c(2000,1500))
# # Make raster stack with two layers according to the environmental variables of the dataframe
# bio <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=2),xmn=0,xmx=2,ymn=0,ymx=1),
#                   "BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=2),xmn=0,xmx=2,ymn=0,ymx=1)))

# Or load ones (Camargue)
load("bio.RData")

###### Genetic data :

# load fake GeneticData 
genetic <- read.table("GeneticData.txt", header=TRUE)
genetic[,"x"] <- rep(0.5, times = 10)
genetic[,"y"] <- rep(0.5, times = 10)

# genetic value of the ancestor
ancestor <- 200

# assuming we have the step values for each locus
steps <- c(1,2,3,4,5,4,4,4,4,4)
length(steps)

# ParamList <- askListOfParameters(rasterStack=bio, nb_simulations=100)
# save(ParamList, file = "Exemples/OnePop/ParamList.RData")

# Or load it from working directory
load("ParamList.RData")

# launch simulations
simSpatialCoal(nbSimul=10000, ParamList=ParamList, rasterStack=bio, GeneticData=genetic, initialGenetValue=ancestor,
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
