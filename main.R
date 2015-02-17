### Setting working directories
rm(list=ls())
# wd="/media/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # portable
# wd="/media/dupas/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # fixe
# wd="/home/legs/GraphPOP/" # portable steph
# wd="/home/dupas/GraphPOP/" # fixe
# wd="/home/arno/These/GraphPOP" # portable arno
# wd="/home/arnaudb/Documents/GraphPOP" # labo arno
# setwd(wd)

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
library(parallel)

########### Parameters initialisation  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

###### Environmental data of temperature and precipitations
# Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900)) 
Data2 <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 
# Make raster stack with two layers according to the environmental variables of the dataframe
rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1),
                          "BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)))

###### Genetic data :

# load fake GeneticData 
GeneticData <- read.table("GeneticData.txt", header=TRUE)

# genetic value of the ancestor
initialGenetValue <- 200

# number of loci under study:
locusNames <- colnames(GeneticData)[!(colnames(GeneticData)%in%c("x","y","Cell_numbers"))]
numberOfLoci <- length(locusNames)

# assuming we have the step values for each locus
stepValueOfLoci <- c(1,2,3,4,5,4,4,4,4,4,2,2,2,2)

# where are the sampled data ?
localizationData <- cellFromXY(object = rasterStack, xy = GeneticData[, c("x", "y")])
names(localizationData)=1:length(localizationData)


###### Asking List to the user

# ParamList <- askListOfParameters(rasterStack=rasterStack, nb_simulations=10000)
# save(ParamList, file = "ParamList.RData")

# Or load it from working directory
load("ParamList.RData")

########## end of parameters initialisation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

local({
  
  # open a connection to a temporary file : pipe between master process and child
  f <- fifo(tempfile(), open="w+b", blocking=T)
  
  if (inherits(parallel:::mcfork(), "masterProcess")) {
    # Child
    progress <- 0.0
    
    while (progress < 1 && !isIncomplete(f)) {
      msg <- readBin(f, "double")
      progress <- progress + as.numeric(msg)
      # send a message in C-style
      cat(sprintf("Progress: %.2f%%\n", progress * 100))
    } 
    
    # close the current child process, informing master process
    parallel:::mcexit()
  }
    
  numJobs <- 1000
  
  mclapply(X = 1:numJobs, FUN = function(x, 
                                                   ParamList, 
                                                   rasterStack, 
                                                   GeneticData, 
                                                   initialGenetValue, 
                                                   numberOfLoci, 
                                                   stepValueOfLoci,
                                                   localizationData){
    
    geneticResults <- matrix(data=NA, nrow=nrow(GeneticData), ncol=numberOfLoci)
    
    ### LOOP ON LOCI >>>>>>>>>>>>>>>>>
    for(locus in 1:numberOfLoci){ # locus=1
      
      # Get the stepValue of the locus under concern
      stepValue <- stepValueOfLoci[locus]
      
      # Get the carrying capacity map :
      rasK <- nicheFunctionForRasterStack(functionList = getFunctionListNiche(ParamList = ParamList, sublist="NicheK"), 
                                          rasterStack = rasterStack,
                                          args = getArgsListNiche(simulation = x, ParamList = ParamList, sublist="NicheK"))
      
      # Get growth rate map :
      rasR <- nicheFunctionForRasterStack(functionList = getFunctionListNiche(ParamList = ParamList, sublist="NicheR"), 
                                          rasterStack = rasterStack,
                                          args = getArgsListNiche(simulation = x, ParamList = ParamList, sublist="NicheR"))
      
      # Get migration matrix :
      kernelMatrix <- dispersionFunctionForRasterLayer(dispersionFunction=getFunctionDispersion(ParamList),
                                                       rasterLayer=rasterStack[[1]], 
                                                       args=getArgsListDispersion(simulation = x, ParamList = ParamList))
      
      migrationMatrix <- migrationRateMatrix(kernelMatrix)
      
      # Get transition matrix :
      transitionBackward <- transitionMatrixBackward(r = values(rasR), K = values(rasK), migration = migrationMatrix)
      transitionForward <- transitionMatrixForward(r = values(rasR), K = values(rasK), migration = migrationMatrix, meth = "non_overlap")
      
      # launch the coalescent
      Coalescent_genetics <- simul_coalescent_only(tipDemes = localizationData,
                                                   transitionForward = transitionForward, 
                                                   transitionBackward = transitionBackward, 
                                                   K = values(rasK))
      
      # adding branch length and genetic data
      Coalescent_genetics <- add_br_length_and_mutation(coalescent = Coalescent_genetics, 
                                                        mutation_rate = ParamList[["Mutation"]][["mutationRate"]][["Values"]][x])
      
      # Transforming the coalescent list into a table
      coalTable <- coalist_2_coaltable(Coalescent_genetics[[1]])
      
      # add resultant 
      whichzero <- coalTable[["mutations"]]==0
      coalTable[["Resultant"]] <- 0
      coalTable[["Resultant"]][!whichzero] <- resultantFunction(nbrMutations = coalTable[["mutations"]][!whichzero],
                                                                stepValue = stepValue,
                                                                mutationModel = getFunctionMutation(ParamList = ParamList),
                                                                args = getArgsListMutation(simulation = x, ParamList = ParamList ))
      
      # add genetic values
      coalTable <- addGeneticValueToCoaltable(coalTable = coalTable, initialGenetValue = initialGenetValue, stepValue = stepValue)
      
      # Record the genetic data
      geneticResults[,locus] <- coalTable[coalTable$coalescing%in%names(localizationData),"genetic_value"]
      # Record the forward log probability
      
    } # END OF LOOP OVER LOCI <<<<<<<<<<<<<
    
    # write results of genetic data 
    fname = paste(getwd(),"/SimulResults/", "Genetics_", x , ".txt", sep="")
    write.table(geneticResults, file=fname)
    write(Coalescent_genetics$forward_log_prob,file=fname,append=TRUE)
    
    # Send progress update
    writeBin(1/numJobs, f)
    
  } # END OF FUNCTION IN MCLAPPLY
  ,
  ParamList = ParamList, 
  rasterStack = rasterStack, 
  GeneticData = GeneticData, 
  initialGenetValue = initialGenetValue, 
  numberOfLoci = numberOfLoci,
  stepValueOfLoci = stepValueOfLoci,
  localizationData = localizationData, mc.cores = 40)
  
  close(f)
  
})

cat("Done\n")

