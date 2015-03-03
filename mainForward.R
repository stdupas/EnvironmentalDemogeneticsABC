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
Dates <- as.Date(startingDate:stoppingDate,origin="1970-01-01");length(Dates)
EnvData = array(c(rep( sample((12:57)*10,400,replace=TRUE),731),rep(sample((15:290)*10,400,replace=TRUE),731)),dim=c(400,731,2),dimnames = list(1:400,as.character(Dates),c("BIO1","BIO12")))
rasterStack <- stack(list("BIO1"=raster(matrix(EnvData[,1,1],nrow=20,ncol=20),xmn=0,xmx=20,ymn=0,ymx=20),
                          "BIO12"=raster(matrix(EnvData[,1,2],nrow=20,ncol=20),xmn=0,xmx=20,ymn=0,ymx=20)))
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

# make Daily raster stack


###### release data 

release <- data.frame(individualNb=1:10,
                      x=runif(10,min=bbox(rasterStack)[1,1],
                              max=bbox(rasterStack)[1,2]),
                      y=runif(10,min=bbox(rasterStack)[2,1],
                              max=bbox(rasterStack)[2,2]),
                      birthDate=as.Date(as.Date("2001/01/03"):as.Date("2001/01/12"),origin="1970-01-01")
                      )
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

fwParamList$NicheR$BIO1$slope$max<-0.03
fwParamList$NicheR$BIO1$slope$Values=runif(10,0,.03)
fwParamList$NicheR$BIO12$slope$max<-0.003
fwParamList$NicheR$BIO12$slope$Values=runif(10,0,.003)

fwParamList$NicheK$BIO1$slope$max<-.1
fwParamList$NicheK$BIO1$slope$Values=runif(10,0,.1)
fwParamList$NicheK$BIO12$slope$max<-0.03
fwParamList$NicheK$BIO12$slope$Values=runif(10,0,.03)
save(list=fwParamList,file="fwParamList.RData")

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
    
  numJobs <- 4
  
  mclapply(X = 1:numJobs, FUN = function(x, 
                                         fwParamList, 
                                         rasterStack, 
                                         individuals, 
                                         localizationData,
                                         EnvData){
    kernelMatrix <- dispersionFunctionForRasterLayer(dispersionFunction=getFunctionDispersion(fwParamList),
                                                     rasterLayer=rasterStack[[1]], 
                                                     args=getArgsListDispersion(simulation = 1, ParamList = fwParamList))
    migrationMatrix <- forwardMigrationRateMatrixFromKernel(kernelMatrix)
#    rK <- generate_parameterSeries(EnvData,
#                                   migrationMatrix,
#                                   nicheKModelsList=getFunctionListNiche(ParamList = fwParamList, sublist="NicheK"),
#                                   nicheRModelsList=getFunctionListNiche(ParamList = fwParamList, sublist="NicheR"),
#                                   paramKList=getArgsListNiche(simulation = x, ParamList = fwParamList, sublist="NicheK"),
#                                   paramRList=getArgsListNiche(simulation = x, ParamList = fwParamList, sublist="NicheR"))
      
    # LOOP OVER DAYS
    for (Date in dimnames(EnvData)[[2]]) #Date=dimnames(EnvData)[[2]][3]
    {
      values(rasterStack) <- EnvData[,Date,]
      # Get the carrying capacity map of the day:
      rasK <- nicheFunctionForRasterStack(functionList = getFunctionListNiche(ParamList = fwParamList, sublist="NicheK"), 
                                          rasterStack = rasterStack,
                                          args = getArgsListNiche(simulation = x, ParamList = fwParamList, sublist="NicheK"))
      indivOfDate <- which(individuals$birthDate==Date)
      # individual in demes exceeding carrying capacity at that date are sampled and removed from the table
      toRemove = individualsToRemoveFromCompet(individuals=indivOfDate,
                                     demeNb=individuals[indivOfDate,"demeNb"],
                                     K=values(rasK))
      if (length(toRemove)>0) individuals = individuals[-toRemove,]
      dateSubset = subset(individuals,individuals$birthDate==Date)
      # launch the siulation
      
      # Get growth rate map :
      rasR <- nicheFunctionForRasterStack(functionList = getFunctionListNiche(ParamList = fwParamList, sublist="NicheR"), 
                                          rasterStack = rasterStack,
                                          args = getArgsListNiche(simulation = x, ParamList = fwParamList, sublist="NicheR"))
      
      rownames(individuals) <- NULL
      newIndividuals <- reprMigr(dateSubset,
               migrationMatrix,
               r=values(rasR),
               generationTime=fwParamList$generationTime$mean$values[x],
               generationTimeRelativeSD=fwParamList$generationTime$SD$values[x])
      if (nrow(newIndividuals)!=0) {rownames(newIndividuals) <- (nrow(individuals)+1):(nrow(individuals)+nrow(newIndividuals))}
      individuals <- rbind(individuals,newIndividuals)
    } # END OF LOOP OVER DAYS <<<<<<<<<<<<<
    
    # write results of genetic data 
    fname = paste(getwd(),"/ForwardSimulResults/", "Individuals_", x , ".txt", sep="")
    write.table(individuals, file=fname)
    
    # Send progress update
    writeBin(1/numJobs, f)
    
  } # END OF FUNCTION IN MCLAPPLY
  ,
  fwParamList = fwParamList, 
  rasterStack = rasterStack, 
  individuals = individuals, 
  localizationData = localisationData, 
  EnvData = EnvData, mc.cores = 2)
  
  close(f)
  
})

cat("Done\n")

