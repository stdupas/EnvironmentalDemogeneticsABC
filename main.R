
### Setting working directories
rm(list=ls())
wd="/media/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # portable
wd="/media/dupas/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # fixe
wd="/home/legs/GraphPOP/" # portable steph
wd="/home/dupas/GraphPOP/" # fixe
wd="/home/arno/These/GraphPOP" # portable arno
wd="/home/arnaudb/Documents/GraphPOP" # labo arno
setwd(wd)

### Sourcing functions files
source("AskModelsFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("CoalescentFunctions.R")
source("PriorFunctions.R")
source("MarkovProcess.R")
source("GeneticDataSimulation.R")
source("SummaryStats.R")
### Sourcing Libraries
library(raster)
library(ape)
library(stringr)
library(lattice)

########### Parameters initialisation  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

###### Environmental data of temperature and precipitations
# Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900)) 
Data2 <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 
# Make raster stack with two layers according to the environmental variables of the dataframe
rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1),
                          "BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)))

# distance method to be applied to genetic data
DistanceMethod = "DeltaMuDistance"

###### Genetic data :

# load fake GeneticData 
load("GeneticData.RData")
save(ParamList,file="ParamList.RData")

# genetic value of the ancestor
initialGenetValue <- 200
# number of loci under study:
locusNames <- colnames(GeneticData)[!(colnames(GeneticData)%in%c("x","y","Cell_numbers"))]
numberOfLoci <- length(locusNames)

# assuming we have the step values for each locus
stepValueOfLoci <- c(1,2,3,4,5)
# where are the sampled data ?
localizationData <- cellFromXY(object = rasterStack, xy = GeneticData[, c("x", "y")])
names(localizationData)=1:length(localizationData)


###### Asking List to the user
# ParamList <- askListOfParameters(rasterStack=rasterStack, nb_simulations=10)
# save(ParamList, file = "ParamList.RData")

# Or load it from working directory
load("ParamList.RData")

# Simulations number
Simulated = 1:2
# forward log probability vector

########## end of parameters initialisation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



######################### Using Functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#### LOOP ON SIMULATIONS >>>>>>>>>>>>>>>>>>>>>>
for(simulation in Simulated){ # simulation <- 1
  geneticResults <- matrix(data=NA, nrow=nrow(GeneticData), ncol=numberOfLoci)
  
  ### LOOP ON LOCI >>>>>>>>>>>>>>>>>
  for(locus in 1:numberOfLoci){ # locus=1
    
    # Get the stepValue of the locus under concern
    stepValue <- stepValueOfLoci[locus]
    
    # Get the carrying capacity map :
    K <- nicheFunctionForRasterStack(functionList = getFunctionListNiche(ParamList = ParamList, sublist="NicheK"), 
                                     rasterStack = rasterStack,
                                     args = getArgsListNiche(simulation=simulation, ParamList = ParamList, sublist="NicheK"))
    
    # Get growth rate map :
    r <- nicheFunctionForRasterStack(functionList = getFunctionListNiche(ParamList = ParamList, sublist="NicheR"), 
                                     rasterStack = rasterStack,
                                     args = getArgsListNiche(simulation=simulation, ParamList = ParamList, sublist="NicheR"))
    
    # Get migration matrix :
    kernelMatrix <- dispersionFunctionForRasterLayer(dispersionFunction=getFunctionDispersion(ParamList),
                                                     rasterLayer=rasterStack[[1]], 
                                                     args=getArgsListDispersion(simulation = simulation, ParamList = ParamList))
    
    migrationMatrix <- migrationRateMatrix(kernelMatrix)
    
    # Get transition matrix :
    transitionBackward <- transitionMatrixBackward(r = values(r), K = values(K), migration = migrationMatrix)
    transitionForward <- transitionMatrixForward(r = values(r), K = values(K), migration = migrationMatrix, meth = "non_overlap")
    
    # launch the coalescent
    Coalescent_genetics <- simul_coalescent_only(tipDemes = localizationData,
                                            transitionForward = transitionForward, 
                                            transitionBackward = transitionBackward, 
                                            K = values(K))
    
    # adding branch length and genetic data
    Coalescent_genetics <- add_br_length_and_mutation(coalescent = Coalescent_genetics, 
                                                      mutation_rate = ParamList[["Mutation"]][["mutationRate"]][["Values"]][simulation])
    
    # Transforming the coalescent list into a table
    coalTable <- coalist_2_coaltable(Coalescent_genetics[[1]])
    
    # add resultant 
    whichzero <- coalTable[["mutations"]]==0
    coalTable[["Resultant"]] <- 0
    coalTable[["Resultant"]][!whichzero] <- resultantFunction(nbrMutations = coalTable[["mutations"]][!whichzero],
                                                  stepValue = stepValue,
                                                  mutationModel = getFunctionMutation(ParamList = ParamList),
                                                  args = getArgsListMutation(simulation = simulation, ParamList = ParamList ))
    # add genetic values
    coalTable <- addGeneticValueToCoaltable(coalTable = coalTable, initialGenetValue = initialGenetValue, stepValue = stepValue)
    
    # Record the genetic data
    geneticResults[,locus] <- coalTable[coalTable$coalescing%in%names(localizationData),"genetic_value"]
    # Record the forward log probability
    
  } # END OF LOOP OVER LOCI <<<<<<<<<<<<<
  
  # write results of genetic data 
  fname = paste(getwd(),"/SimulResults/", "Genetics_",simulation, ".txt", sep="")
  write.table(geneticResults, file=fname)
  write(coalescentList$forward_log_prob,file=fname,append=TRUE)
} # END OF LOOP OVER SIMULATIONS <<<<<<<<<<<<<<<<<<<

# Post traitement


# Caracterizing rotation to apply to genetic simulations to get summary statistics
#

Rotation = PCA_rotation(GeneticData[,locusNames],DistanceMethod)
summaryStatObsMat = as.matrix(do.call(what = DistanceMethod, args = list(GeneticData[,locusNames])))%*%Rotation
colN <- paste(rep(colnames(summaryStatObsMat),each=length(rownames(summaryStatObsMat))),rownames(summaryStatObsMat),sep=".")
summaryStatObsVect = c(summaryStatObsMat)
names(summaryStatObsVect)<- colN
summaryStatObs <- cbind(ForWLogProb=0,as.data.frame(t(summaryStatObsVect)))
summaryStatSim=summaryStatObs[-1,]
#summaryStatSim <- summaryStatSim[-1,]

## filling summary statistics


Simulated <- grep (pattern = "^Genetics_\\d*.txt$", x=list.files("SimulResults/"), value = TRUE)
for (filen in Simulated)
{
  genetic_i <- read.table(paste("SimulResults/",filen,sep=""),nrows=nrow(GeneticData))
  ForwLogL <- as.numeric(scan(paste("SimulResults/",filen,sep=""),
                   what="numeric",
                   skip=nrow(GeneticData)+1))
  genetDist = do.call(what = DistanceMethod, 
                      args = list(genetic_i))
  rotatedGenetDist = as.matrix(genetDist)%*%Rotation
  i <- grep(pattern="\\d", x=filen, value=TRUE)
  summaryStatSim[i,] <- c(ForwLogL,rotatedGenetDist)
}


