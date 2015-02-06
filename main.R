
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

###### Genetic data :

# load fake GeneticData 
load("GeneticData.RData")

# number of loci under study:
numberOfLoci <- ncol(GeneticData[-c(1,2,3)])
# assuming we have the step values for each locus
stepValueOfLoci <- c(1,2,3,4,5)

# where are the sampled data ?
localizationData <- cellFromXY(object = K, xy = GeneticData[, c("x", "y")])
names(localizationData)=1:length(localizationData)


###### Asking List to the user
# ParamList <- askListOfParameters(rasterStack=rasterStack, nb_simulations=10)
# save(ParamList, file = "ParamList.RData")

# Or load it from working directory
load("ParamList.RData")


########## end of parameters initialisation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



######################### Using Functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#### LOOP ON SIMULATIONS >>>>>>>>>>>>>>>>>>>>>>
for(simulation in 1:2 ){ # simulation <- 1
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
    coalescentList <- simul_coalescent_only(tipDemes = localizationData,
                                            transitionForward = transitionForward, 
                                            transitionBackward = transitionBackward, 
                                            K = values(K))
    
    # adding branch length and genetic data
    Coalescent_genetics <- add_br_length_and_mutation(coalescentList, mutation_rate=.1)
    
    # Transforming the coalescent list into a table
    coalTable <- coalist_2_coaltable(Coalescent_genetics[[1]])
    
    # add resultant 
    coalTable[["Resultant"]] <- resultantFunction(nbrMutations = coalTable[["mutations"]],
                                                  stepValue = stepValue,
                                                  mutationModel = getFunctionMutation(ParamList = ParamList),
                                                  args = getArgsListMutation(simulation = simulation, ParamList = ParamList ))
    
    # add genetic values
    coalTable <- addGeneticValueToCoaltable(coalTable,200,stepValue)
    
    # Record the genetic data
    geneticResults[,locus] <- coalTable[coalTable$coalescing%in%names(localizationData),"genetic_value"]
    
  } # END OF LOOP OVER LOCI <<<<<<<<<<<<<
  
  # write results of genetic data 
  write.table(t(geneticResults), file=paste(getwd(),"/SimulResults/", "Genetics_",simulation, ".txt", sep=""))
  
} # END OF LOOP OVER SIMULATIONS <<<<<<<<<<<<<<<<<<<
