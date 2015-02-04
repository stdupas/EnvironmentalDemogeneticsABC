
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


# where are the sampled data ?
localizationData <- c(1,4,2,2,1,1,2,3);
names(localizationData)=1:length(localizationData)

###### Genetic parameters :
mutation_rate=1E-4

###### Asking List to the user
# ParamList <- askListOfParameters(rasterStack=rasterStack, nb_simulations=10)
# save(ParamList, file = "ParamList.RData")

# Or load it from working directory
load("ParamList.RData")

########## end of parameters initialisation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



######################### Using Functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get the reference table
referenceTable <- referenceTableFromList(ParamList)
# just imagine you pass the reference table into an apply function :
parameters <- referenceTable[,1]

#### Imagine we are looping over simulations of coalescent :
simulation <- 1

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
                                           args=getArgsListDispersion(simulation=simulation, ParamList=ParamList))

migrationMatrix <- migrationRateMatrix(kernelMatrix)

# Get transition matrix :
transitionBackward <- transitionMatrixBackward(r=values(r), K=values(K), migration= migrationMatrix)
transitionForward <- transitionMatrixForward(r=values(r), K=values(K), migration= migrationMatrix, meth="non_overlap")

# launch the coalescent
coalescentList <- simul_coalescent_only(tipDemes=localizationData,transitionForward=transitionForward, transitionBackward=transitionBackward, K= values(K))

# adding branch length and genetic data
Coalescent_genetics <- add_br_length_and_mutation(coalescentList,mutation_rate=.1)
coaltable <- coalist_2_coaltable(Coalescent_genetics[[1]])


