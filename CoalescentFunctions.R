
simul_coalescent <- function(geneticData, rasterStack, pK, pr, shapesK, shapesr, shapeDisp, pDisp,
                             mutation_rate, initial_genetic_value, mutation_model, stepvalue, mut_param)
{
  # Simulates a coalescent in a lansdcape characterized by an environmental variable rasterStack, for a species with a given niche function. 
  #
  # Args:
  #   geneticData: genetic data table : with coordinates
  #   rasterStack: environmental variables raster stack
  #   pK: parameters values of K (Xmin, Xmax, Xopt, YXmin, Yxmax, Yopt) for each environmental variable
  #   pr: parameters values of r (Xmin, Xmax, Xopt, YXmin, Yxmax, Yopt) for each environmental variables
  #   shapesK: shapes of niche function (reaction norm) for the carrying capacity K
  #   shapesr: shapes of niche function (reaction norm) for the growth rate r
  #   pDisp: parameters of dispersion
  #   mutation_rate: mutation rate
  #   initial_genetic_value: genetic value attributed to the common ancestor
  #   mutation_model: specify a model of mutation : "step_wise" (Stepwise Mutation Model),"tmp" (Two Phases Mutation Model)
  #   stepvalue:  
  #   mut_param: 
  # Returns :
  #   A list with all coalescence informations :
  #   List of 4
  #    $ coalescent      :List of "i" (with "i" the number of coalescence events, coalescence involving multiple individuals counts for 1 event.)
  #      ..$ :List of 5
  #      ..  ..$ time           : time of coalescence
  #      ..  ..$ coalescing     : coalescing nodes
  #      ..  ..$ new_node       : "new" in a backward sense, ie the node resulting of the coalescence of the coalescing nodes
  #      ..  ..$ br_length      : the length of the branches
  #      ..  ..$ mutations      : the number of mutations which occured along the branch
  #    $ mutation_rate          : 
  #    $ forward_log_prob       : the logarithm of the forward probability of this coalescent
  #    $ genetic_values         : a data frame with : time, coalescing, new_node, br_length, mutations (0 or NA), genetic_value (initial), resultant (0 or NA)
  
  ### Computes Niche Function variables :
  K = ReactNorm(values(rasterStack),pK,shapesK)[,"Y"]
  r = ReactNorm(values(rasterStack),pr,shapesr)[,"Y"]
  
  ### Computes stochastic matrix :
  migrationM <- migrationMatrix(rasterStack,shapeDisp, pDisp)
  transitionmatrice <- transitionMatrixBackward(r, K, migration= migrationM)
  transition_forward <- transitionMatrixForward(r, K, migration= migrationM)
  
  #### Initialize variables needed for the coalescent simulation process :
  time=0
  prob_forward=NA
  number_of_nodes_over_generations=0
  N <- round(K)
  coalescent = list() 
  # Nodes are initialized : 1 individual <=> 1 node
  nodes = as.numeric(rownames(geneticData))
  names(nodes)=as.character(nodes)
  # Cells in which were the genes sampled in the landscape :
  cell_number_of_nodes <- geneticData[,"Cell_numbers"]
  names(cell_number_of_nodes) <- nodes
  # Cells in which the previous generation genes were in the landscape :
  parent_cell_number_of_nodes <- cell_number_of_nodes 
  # A list of cells with all the genes remaining in each cell. Initialized and optimized.
  nodes_remaining_by_cell = list()
  cell <- as.array(seq(from=1,to=ncell(rasterStack),by=1))
  nodes_remaining_by_cell <- lapply(X=cell, FUN=remainingNodes, cell_number_of_nodes)
  
  # Number of coalescence events :
  single_coalescence_events=0 
  single_and_multiple_coalescence_events=0
  
  ### Simulating the coalescent process :
  while (length(unlist(nodes_remaining_by_cell))>1) 
  {
    
    ## Migration
    # we localize the parents in the landscape by sampling in the backward transition matrix. Optimized.
    names_node <- names(parent_cell_number_of_nodes)
    parent_cell_number_of_nodes <- apply(X=as.array(1:length(parent_cell_number_of_nodes)), MARGIN=1, FUN=backwardParentsLocalizationSampling, rasterStack, transitionmatrice, cell_number_of_nodes) 
    names(parent_cell_number_of_nodes) <- names_node
    # once we know the parent cell numbers, we calculate the forward dispersion probability of the event
    prob_forward[time] = sum(log(transition_forward[parent_cell_number_of_nodes,cell_number_of_nodes]))
    number_of_nodes_over_generations = number_of_nodes_over_generations + length(cell_number_of_nodes)
    
    ## Coalescence
    time=time+1; if (round(time/10)*10==time) {print(time)}
    
    # we now perform coalescence within each cell of the landscape for the parents
    for (cell in 1:ncell(rasterStack))#cell=1;cell=2;cell=3;cell=4;cell=5;cell=26;cell=10
    {
      # add a local variable, easier to manipulate than the reference to a list...
      nodes_remaining_in_the_cell = nodes_remaining_by_cell[[cell]] <- as.numeric(names(which(parent_cell_number_of_nodes==cell)))
      
      # we obtain the identities in the geneticData table (line) of the nodes remaining in the cell
      if (length(nodes_remaining_in_the_cell)>1)
      {
        # Create a function for parentality attribution within a cell :
        parentoffspringmatrix <- parentalityAttributationWithinACell(nodes_remaining_in_the_cell=nodes_remaining_in_the_cell, N=N, cell=cell)
        
        # Columns of parentoffspringmatrix with more than one TRUE allow to identify coalescing individuals :
        if (any(colSums(parentoffspringmatrix)>1) )
        {
          #  Loop over all the parents in which coalescence event occur
          for (multiple in which(colSums(parentoffspringmatrix)>1)) # multiple<-which(colSums(parentoffspringmatrix)>1)[1]
          {
            # Record the coalescence event : 
            single_coalescence_events = single_coalescence_events +1
            
            # which(parentoffspringmatrix[,multiple]) identifies which node in the column coalesce
            nodes_that_coalesce = names(which(parentoffspringmatrix[,multiple]))
            
            # attibutes new node number to the ancestor
            new_node <- max(nodes)+1
            # removes the nodes that coalesced from the node vector
            nodes=nodes[!(names(nodes)%in%nodes_that_coalesce)]
            # adds them to the nodes vector
            nodes=append(nodes,new_node)
            names(nodes)[length(nodes)]=new_node
            
            # updating of vector parent_cell_number_of_nodes (adding the cell number of the new node and removing the nodes that disapeared)
            parent_cell_number_of_nodes <- append(parent_cell_number_of_nodes[!(names(parent_cell_number_of_nodes)%in%nodes_that_coalesce)],cell)
            names(parent_cell_number_of_nodes)[length(parent_cell_number_of_nodes)]<-new_node
            # adds the event to the list coalescent: time, which node coalesced, and the number of the new node
            coalescent[[single_coalescence_events]] <- list(time=time,coalescing=as.numeric(nodes_that_coalesce),new_node=new_node)
            # updating the nodes vector for the cell
            nodes_remaining_in_the_cell = nodes_remaining_by_cell[[cell]] <- append(nodes_remaining_in_the_cell[!nodes_remaining_in_the_cell %in% nodes_that_coalesce],new_node)
            # updates the number of coalescent events 
            single_and_multiple_coalescence_events = single_and_multiple_coalescence_events + length(nodes_that_coalesce) - 1
            
          } #  end of loop over all the parents in which coalescence event occur
        } # end of the if condition "there are coalescing events"
      } # end of the condition "there are more than 1 individual in the cell
    } # end of the loop across the cells
    
    cell_number_of_nodes = parent_cell_number_of_nodes
  } # end of the backward generation while coalescence loop
  
  # Inform historical and mutational processes
  coalescent=add_br_length_and_mutation(coalescent,mutation_rate,initial_genetic_value)
  # Formatting the output
  list(coalescent=coalescent,mutation_rate=mutation_rate,
       forward_log_prob=sum(prob_forward)/number_of_nodes_over_generations,
       genetic_values=genetics_of_coaltable(coalist_2_coaltable(coalescent),initial_genetic_value,mutation_model,stepvalue,mut_param))
  # forward_log_prob is the average per generation of the log probability of the forward movements of the genes in the landscape
}

parentalityAttributationWithinACell <- function(nodesRemainingInCell, N, cell)
{
  # Create a function for parentality attribution within a cell used in simul_coalescent
  #
  # Args:
  #   nodesRemainingInCell:
  #   N: a vector giving the population size in each cell
  #   cell: numeric giving the cell we want to analyse
  #
  # Returns:
  #   A matrix of parentality attribution. Two nodes coalesce if they have TRUE for the same parent (parents are in columns)
  
  nbGenesRemaining=length(nodesRemainingInCell)
  
  # Attribute parents (among K possible parents) to each node present in the cell
  smp = sample(N[cell],length(nodesRemainingInCell),replace=TRUE)
  
  # A logical matrix in which lines represent the nodes in the cell and column represent their parent :
  # (actually, this line is a simple test to transform the parentality info under a TRUE/FALSE form)
  parentoffspringmatrix <- matrix(smp, nrow=nbGenesRemaining, ncol=N[cell]) == matrix(1:N[cell], nrow=nbGenesRemaining, ncol=N[cell], byrow=TRUE)
  rownames(parentoffspringmatrix) <- nodesRemainingInCell
  
  return(parentoffspringmatrix)
}

backwardParentsLocalizationSampling <- function(node, rasterStack, transitionMatrice, CellIdOfNodes)
{
  # Localizes the parents in the landscape by sampling in the backward transition matrix.
  #
  # Args:
  #   node: the node from which to set parents localization
  #   rasterStack: the rasterStack used to define landscape structure (cells)
  #   transitionMatrice: the transition matrice giving the probability the migrate between cells of the landscape
  #   CellIdOfNodes: the vector giving nodes localization (giving the number of the cell where each node sits)
  #
  # Returns:
  #   The number of the cell where the parents are likely to be found
  return(sample(ncell(rasterStack),size=1,prob=c(transitionMatrice[CellIdOfNodes[node],])))
}

remainingNodes <- function(cell, cellIdOfNodes)
{
  # Finds all the nodes remaining a specified cell
  #
  # Args:
  #   cell: a numerci giving the cell in which we want to find the remaining nodes
  #   cellIdOfNodes: gives the nodes localization
  return(which(cellIdOfNodes==cell))
}

coalist_2_coaltable <- function(coalist)
{
  coaldf <- data.frame(Reduce(rbind,coalist))
  coaltable <- coaldf[rep(1:dim(coaldf)[1],unlist(lapply(coaldf$coalescing, length))),]
  coaltable[,c("coalescing","br_length","mutations")] <- unlist(coaldf[,c("coalescing","br_length","mutations")])
  coaltable
}
