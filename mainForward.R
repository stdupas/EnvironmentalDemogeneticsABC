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
#source("MutationFunctions.R")
#source("CoalescentFunctions.R")
source("PriorFunctions.R")
#source("MarkovProcess.R")
#source("GeneticDataSimulation.R")
source("forwardModels.R")

### Sourcing Libraries
library(raster)
library(ape)
library(stringr)
library(lattice)
library(parallel)
#library(lubridate)

########### Parameters initialisation  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

###### Environmental data of temperature and precipitations
# Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900)) 
#Data2 <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 
startingDate = as.Date("2001-01-01")
stoppingDate = as.Date("2003-01-01")
Dates <- as.Date(startingDate:stoppingDate);length(Dates)
EnvData = array(c(rep(c(300,120,120,400),731),rep(c(2000,350,350,2900),731)),dim=c(1,4,731,2),dimnames = list(1,1:4,as.character(Dates),c("BIO1","BIO12")))
variables=c("BIO1","BIO12")
# Make raster stack with two layers according to the environmental variables of the dataframe
#rS=stack(raster(as.matrix(EnvData[,,1,1])))
#for (variable in variables)
#{
#  for (Date in as.character(Dates))
#  {
#    rS <- addLayer(rS, raster( as.matrix(EnvData[,,Date,variable] ) ))
#  }
#}

rasterStack <- stack(list("BIO1"=raster(as.matrix(t(EnvData[,,1,1])),xmn=0,xmx=4,ymn=0,ymx=1),
                          "BIO12"=raster(as.matrix(t(EnvData[,,1,2])),xmn=0,xmx=4,ymn=0,ymx=1)))
# make Daily raster stack


###### release data 

release <- data.frame(individualNb=1:10,
                      x=runif(10,min=bbox(rasterStack)[1,1],
                              max=bbox(rasterStack)[1,2]),
                      y=runif(10,min=bbox(rasterStack)[2,1],
                              max=bbox(rasterStack)[2,2]),
                      birthDate=as.Date("2013/01/23"))

individuals = release[,c("individualNb","birthDate")]
individuals$demeNb <- cellFromXY(object = rasterStack, xy = release[, c("x", "y")])

# where are the sampled data ?
localizationData <- individuals$demeNb
names(localizationData)=individuals$individualNb


###### Asking List to the user

# fwParamList <- askListOfParameters(rasterStack=rasterStack, nb_simulations=10)
# fwParamList$generationTime$mean$priorLaw <- "uniform"
# ParamList$generationTime$mean$min <- 15
# ParamList$generationTime$mean$max <- 30
# ParamList$generationTime$mean$values <- runif(10,15,30)
# ParamList$generationTime$SD$priorLaw <- "uniform"
# ParamList$generationTime$SD$min <- .1
# ParamList$generationTime$SD$max <- .1
# ParamList$generationTime$SD$values <- rep(.1,10)
#


# save(fwParamList, file = "fwParamList.RData")

# Or load it from working directory
load("fwParamList.RData")

########## end of parameters initialisation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# simulation 1
x=1


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
    
  numJobs <- 10
  
  mclapply(X = 1:numJobs, FUN = function(x, 
                                         fwParamList, 
                                         rasterStack, 
                                         individuals, 
                                         localizationData,
                                         EnvData){
    # LOOP OVER DAYS
    for (Date in dimnames(EnvData)[[3]]) Date=dimnames(EnvData)[[3]][1]
    {
      
      # launch the siulation
      # Get the carrying capacity map :
      rasK <- nicheFunctionForArray(functionList = getFunctionListNiche(ParamList = fwParamList, sublist="NicheK"), 
                                    array = EnvData[,,Date,],
                                    args = getArgsListNiche(simulation = x, ParamList = fwParamList, sublist="NicheK"))
      
      # Get growth rate map :
      rasR <- nicheFunctionForRasterStack(functionList = getFunctionListNiche(ParamList = fwParamList, sublist="NicheR"), 
                                          rasterStack = rasterStack,
                                          args = getArgsListNiche(simulation = x, ParamList = fwParamList, sublist="NicheR"))
      
      # Get migration matrix :
      kernelMatrix <- dispersionFunctionForRasterLayer(dispersionFunction=getFunctionDispersion(ParamList),
                                                       rasterLayer=rasterStack[[1]], 
                                                       args=getArgsListDispersion(simulation = x, ParamList = fwParamList))
      migrationMatrix <- forwardMigrationRateMatrixFromKernel(kernelMatrix)
      
      # Get transition matrix :
      absTransForw <- absoluteForwardTransition(r = values(rasR), K = values(rasK), migration = migrationMatrix)
      
    } # END OF LOOP OVER DAYS <<<<<<<<<<<<<
    
    # write results of genetic data 
    fname = paste(getwd(),"/SimulResults/", "Genetics_", x , ".txt", sep="")
    write.table(geneticResults, file=fname)
    write(Coalescent_genetics$forward_log_prob,file=fname,append=TRUE)
    
    # Send progress update
    writeBin(1/numJobs, f)
    
  } # END OF FUNCTION IN MCLAPPLY
  ,
  fwParamList = fwParamList, 
  rasterStack = rasterStack, 
  GeneticData = GeneticData, 
  initialGenetValue = initialGenetValue, 
  numberOfLoci = numberOfLoci,
  stepValueOfLoci = stepValueOfLoci,
  localizationData = localizationData, mc.cores = 2)
  
  close(f)
  
})

cat("Done\n")

