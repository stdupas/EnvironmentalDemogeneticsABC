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
  #   the files of simulated genetic values in the SimulResults repertory
  
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

      
      ### LOOP ON LOCI >>>>>>>>>>>>>>>>>
      
      for(locus in 1:numberOfLoci){ # locus=1
        
        # Get the stepValue of the locus under concern
        stepValue <- stepValueOfLoci[locus]
        
        maxCoalEvent <- length(localizationData) - 1
        
        # coalescent informations : (time of coalescence, Child 1, Child 2, Parent, Branch Length, mutation nbr, resultant, genet values)
        coal <- matrix(data = NA, nrow = maxCoalEvent, ncol = 8)    
        
        # launch the coalescent
        coal[,c(1:4)] <- spatialCoalescentSimulation(tipDemes = localizationData, 
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


spatialCoalescentSimulation <- function(tipDemes, transitionBackward, N){
  # Simulate a genealogy backward in the time, accross demes
  # 
  # Args:
  #   tipdDemes: vector of the demes in which each node is found a time 0.
  #   transitionBackward: matrix of transition backward in time
  #   N a vector of population sizes
  #
  # Returns: 
  #   A matrix describing the coalescence events : time/childNode1/childNode2/parentNode
  
  ###### INITIALISATION
  time <- 0
  events <- 0
  headNode <- length(tipDemes)
  maxCoalEvent <- length(tipDemes) - 1
  nodesState <- c(tipDemes, rep(NA, maxCoalEvent))
  
  # coalescent informations : (time of coalescence, Child 1, Child 2, Parent)
  coalescent <- matrix(data = NA, nrow = maxCoalEvent, ncol = 4)
  
  ###### REPEAT UNTIL TOTAL COALESCENCE
  while (is.na(tail(nodesState, n=1))){
    time <- time +1
    
    #### MIGRATION
    nodesState[!is.na(nodesState)] <- vapply(X = nodesState[!is.na(nodesState)],
                                             FUN = function(x, N, transitionBackward)
                                             {sample( length(N), size = 1, prob = c(transitionBackward[x,]) )},
                                             N = N, transitionBackward = transitionBackward,
                                             FUN.VALUE = c(1))
        
    ####### CANDIDATES NODES FOR COALESCENCE
    # for active nodes, i.e which are not coded by NA : 
    activeNodes <- which(!is.na(nodesState))
    activeDemes <- nodesState[activeNodes]
    # gives indices of the demes that are duplicated
    dup <- which(duplicated(activeDemes) | duplicated(activeDemes, fromLast= TRUE))
    # gives the demes in which more than one node exist :
    demes <- unique(activeDemes[dup])
    # gives a list of nodes who can perhaps coalesce 
    candidates <- lapply(X = demes,
                         FUN = function(x, nodesState){which(nodesState == x)},
                         nodesState = nodesState)
    
    ####### COALESCENCE
    if(length(candidates) > 0){
      
      for(x in seq(from = 1, to = length(candidates))){ # x <- 1
        
        focalDeme <- demes[x]
        # /!\ If N=0, the nodes would automatically coalesce (parents n°0 for everyone) -> make sure this does not happen !
        if(N[focalDeme]==0){stop(paste("in spatialCoalescentSimulation you are trying to coalesce in an empty deme : in deme", x ,", N=0"))}
        
        # Attribute parents (among N possible parents) to each node present in the deme
        parents <- sample(N[focalDeme], size = length(candidates[[x]]), replace = TRUE) # parents[1] <- parents[2]
        # Test for equality of parents :
        anonymous <- which(duplicated(parents) | duplicated(parents, fromLast= TRUE))
        
        # If nodes have same parent node
        if(length(anonymous) > 1) {
          
          # sample in the candidates nodes who will coalesce
          children <- sample(x = candidates[[x]], size = length(anonymous), replace = FALSE)
          # number of new coalescent events
          nEvents <- length(children) -1
          
          # Move header node forward, and skip ephemeral ones
          headNode <- headNode + nEvents
          # Precise the deme were araised the new node
          nodesState[headNode] <- focalDeme
          # Shut down children nodes
          nodesState[children] <- NA
          
          lines <- seq(from = events+1, to = events + nEvents)
          parentNodes <- seq(from = headNode - nEvents + 1, to = headNode )
          # Fill time
          coalescent[lines, 1] <- rep(x = time, times = nEvents)
          # Fill Child1
          coalescent[lines, 2] <- c(children[1], parentNodes[-length(parentNodes)])
          # Fill Child2
          coalescent[lines, 3] <- c(children[-1])
          # Fill parents
          coalescent[lines, 4] <- parentNodes
          events <- events + nEvents
          
        } # end of if there are coaelescing nodes
      } # end of for loop over demes
    } # end of if there are co occuring nodes in the same deme
  } # end of while coalescence is not complete
  return(coalescent)
}


coalescent_2_newick <- function(coalescent)
{
  # coalescent_2_newick
  # function that converts coalescent to newick format tree
  # argument: coalescent list 
  # value : newwick foramt tree
  # 
  # Example
  # trB = matrix(c(1/4,1/2,1/4,0,1/3,1/3,1/6,1/6,1/2,1/4,1/8,1/8,1/5,2/5,2/5,0),nrow=4,ncol=4,byrow=TRUE)
  # trF = matrix(c(1/4,1/2,1/2,1/8,1/3,1/3,1/6,1/6,1/2,1/4,1/8,1/8,0,0,2/5,0),nrow=4,ncol=4,byrow=TRUE)
  # K=c(4,3,1,5)
  # tipsDemes = c(1,4,2,2,1,1,2,3);names(tipsDemes)=1:8
  # Coalescent = simul_coalescent_only(tipDemes=tipsDemes,transitionForward=trF,transitionBackward=trB,K=K)
  # coalescent_2_newick(Coalescent)

  tree=paste(" ",coalescent[[1]][[length(coalescent)]]$new_node," ",sep="")
  for (i in length(coalescent):1)
  {
    Time = coalescent[[1]][[i]]$time
    coalesc <- as.character(coalescent[[1]][[i]]$coalescing)
    tree <- str_replace(tree,paste(" ",as.character(coalescent[[1]][[i]]$new_node)," ",sep=""),paste(" ( ",paste(" ",coalesc," :",coalescent[[1]][[i]]$br_length,collapse=" ,",sep=""),") ",sep=""))
  }
  tree <- gsub(" ","",paste(tree,";",sep=""))
  tree
}


plotCoalescentGenetics <- function(coalescent,genetic_table,with_landscape=FALSE,legend_right_move=-.2)
{
  # function that plots a coalecent, with tips demes as specific color
  # argument: coalescent list 
  #
  # 
  # Example
  # trB = matrix(c(1/4,1/2,1/4,0,1/3,1/3,1/6,1/6,1/2,1/4,1/8,1/8,1/5,2/5,2/5,0),nrow=4,ncol=4,byrow=TRUE)
  # trF = matrix(c(1/4,1/2,1/2,1/8,1/3,1/3,1/6,1/6,1/2,1/4,1/8,1/8,0,0,2/5,0),nrow=4,ncol=4,byrow=TRUE)
  # K=c(4,3,1,5)
  # tipsDemes = c(1,4,2,2,1,1,2,3);names(tipsDemes)=1:8
  # Coalescent = simul_coalescent_only(tipDemes=tipsDemes,transitionForward=trF,transitionBackward=trB,K=K)
  # plotCoalesentGenetics(coalescent_2_newick(Coalescent),tipDemes,legend_right_move=.2)
  
  par(mfrow=c(1,1),oma=c(0,0,0,4),xpd=TRUE)
  tipcells <- tipDemes[as.numeric(read.tree(text=coalescent_2_newick(coalescent))$tip.label)]
   #tipcells <- geneticData$Cell_numbers[as.numeric(coalescent_2_phylog(coalescent)$tip.label)]
  tipcols = rainbow(ncell(rasK))[tipcells]
  phylog_format_tree <- coalescent_2_phylog(coalescent)
  phylog_format_tree$tip.label <- paste(phylog_format_tree$tip.label,genetic_table[order(genetic_table$coalescing)[as.numeric(phylog_format_tree$tip.label)],"genetic_value"],sep=":")
  plot(phylog_format_tree,direction="downward",tip.color=tipcols)
  legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(legend_right_move,0))
  if (with_landscape) {plot(rasK)}
}