simSpatialCoal <- function(nbSimul, ParamList, rasterStack, GeneticData, initialGenetValue, stepValueOfLoci, cores){
  # Estimates the parameters of a model (spatial, niche, coalescence) in an abc framework
  #
  # Args:
  #   nbSimul: the number of simulations wanted for abc estimation. need to be equal to the number of simulation specified in ParamList
  #   ParamList: the R object describing the model, constructed by askListOfParameters function
  #   rasterStack: the raster object describing environment used for niche modelling and dispersion computations
  #   GeneticData: a matrix giving in row the individuals, in columns the coordinates and the loci : names and order have to be : x, y, ... and names of loci
  #   initialGenetValue: a vector giving the genetic value attributed to the ancestor gene.
  #   stepValueOfLoci: a vector giving the assumed step value for each locus, given in the same order as in GeneticData
  #   cores: the number of cores used for computation
  #
  # Returns:
  #   the results of the abc analysis, and the files of simulated genetic values in the SimulResults repertory
   
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
  
  # Create a directory to store simulations results
  dir.create(path=paste(getwd(), "/SimulResults", sep=""))
  
  # number of loci under study:
  locusNames <- colnames(GeneticData)[!(colnames(GeneticData)%in%c("x","y"))]
  numberOfLoci <- length(locusNames)
  
  # where are the sampled data ?
  localizationData <- cellFromXY(object = rasterStack, xy = GeneticData[, c("x", "y")])
  #names(localizationData)=1:length(localizationData)
    
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
    
    numJobs <- nbSimul
    
    ######## Launch the simulations in parallel ############################################### 
    
    mclapply(X = 1:numJobs, FUN = function(x, 
                                           ParamList, 
                                           rasterStack, 
                                           GeneticData, 
                                           initialGenetValue, 
                                           numberOfLoci, 
                                           stepValueOfLoci,
                                           localizationData){
      
      geneticResults <- matrix(data=NA, nrow=nrow(GeneticData), ncol=numberOfLoci)
      forwardProb <- c()
      
      # Get the mutation rate
      mutationRate <- ParamList[["Mutation"]][["mutationRate"]][["Values"]][x]
      
      # Get the carrying capacity map :
      rasK <- nicheFunctionForRasterStack(functionList = getFunctionListNiche(ParamList = ParamList, sublist="NicheK"), 
                                          rasterStack = rasterStack,
                                          args = getArgsListNiche(simulation = x, ParamList = ParamList, sublist="NicheK"))
      rasK <- round(rasK)
      
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
      
      ### LOOP ON LOCI >>>>>>>>>>>>>>>>>
      
      for(locus in 1:numberOfLoci){ # locus=1
        
        # Get the stepValue of the locus under concern
        stepValue <- stepValueOfLoci[locus]
        
        maxCoalEvent <- length(localizationData) - 1
        
        # coalescent informations : (time of coalescence, Child 1, Child 2, Parent, Branch Length, mutation nbr, resultant, genet values)
        coal <- matrix(data = NA, nrow = maxCoalEvent, ncol = 8)    
        
        # launch the coalescent
        coal[,c(1:4)] <- spatialCoalescentSimulation(tipDemes = localizationData, transitionForward = transitionForward, 
                                                           transitionBackward = transitionBackward, 
                                                           N = round(values(rasK)))
        
        # add branch length
        coal[,5] <- c(coal[,1][-1] - coal[,1][-nrow(coal)], NA)
        
        # add mutation number
        coal[,6] <- vapply(X = coal[,5],
                           FUN = function(x){rpois(n = 1, lambda = mutationRate*x)},
                           FUN.VALUE = c(1))
        
        # add resultant 
        coal[,7] <- resultantFunction(nbrMutations = coal[,6],
                                      stepValue = stepValue,
                                      mutationModel = getFunctionMutation(ParamList = ParamList),
                                      args = getArgsListMutation(simulation = x, ParamList = ParamList ))
        
        # add genetic values
        coal[nrow(coal),8] <- initialGenetValue
        for(i in seq(from = nrow(coal)-1, to = 1)){ coal[i,8] <- coal[i+1,8] + coal[i,7] }
        
        # Record the genetic data
        n2 <- which(coal[,2] %in% seq(from = 1, to = length(localizationData)))
        n3 <- which(coal[,3] %in% seq(from = 1, to = length(localizationData)))
        geneticResults[coal[n2, 2], locus] <- coal[n2, 8]
        geneticResults[coal[n3, 3], locus] <- coal[n3, 8]
        
        
      } # END OF LOOP OVER LOCI <<<<<<<<<<<<<
      
      # write results of genetic data 
      fname = paste(getwd(),"/SimulResults/", "Genetics_", x , ".txt", sep="")
      write.table(geneticResults, file=fname)
      
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
    localizationData = localizationData, 
    mc.cores = cores,
    mc.preschedule = FALSE)
    
    close(f)
    
  })
  
  cat("Simulations Done\n")

}



