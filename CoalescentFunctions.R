
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
  cell <- as.array(seq(from=1,to=length(K),by=1))
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
    time=time+1; 
    
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

simul_coalescent_only <- function(tipDemes,transitionForward, transitionBackward, K)
{
  # Simulates a coalescent in a lansdcape characterized by an environmental variable rasterStack, for a species with a given niche function. 
  #
  # Args:
  #   tipDemes : deme number in the lanscape  of the individuals we study the coalescence
  #                       (corresponds to lines and columns in the transition matrixes)
  #   K : vector of carrying capacity for coalescence probability calculation
  #   transitionForward : forward transtion probability matrix among demes to calculate forward likelihood
  #   transitionBackward : backward transition probability matrix among demes to simulate backard movements
  #
  # Returns :
  #   A list with all coalescence informations :
  #   List of 4
  #    $ coalescent      :List of "i" (with "i" the number of coalescence events, coalescence involving multiple individuals counts for 1 event.)
  #      ..$ :List of 5
  #      ..  ..$ time           : time of coalescence
  #      ..  ..$ coalescing     : coalescing nodes
  #      ..  ..$ new_node       : "new" in a backward sense, ie the node resulting of the coalescence of the coalescing nodes
  # Example
  # trB = matrix(c(1/4,1/2,1/4,0,1/3,1/3,1/6,1/6,1/2,1/4,1/8,1/8,1/5,2/5,2/5,0),nrow=4,ncol=4,byrow=TRUE)
  # trF = matrix(c(1/4,1/2,1/2,1/8,1/3,1/3,1/6,1/6,1/2,1/4,1/8,1/8,0,0,2/5,0),nrow=4,ncol=4,byrow=TRUE)
  # K=c(4,3,1,5)
  # tipsDemes = c(1,4,2,2,1,1,2,3);names(tipsDemes)=1:8
  # simul_coalescent_only(tipDemes=tipsDemes,transitionForward=trF,transitionBackward=trB,K=K)

  #### Initialize variables needed for the coalescent simulation process :
  
  time=0
  prob_forward=NA
  number_of_nodes_over_generations=0
  coalescent = list() 
  # Nodes are initialized : 1 individual <=> 1 node
  nodes = as.numeric(names(tipDemes))
  names(nodes)=as.character(nodes)
  # initiation of vectors hosting the deme number of descendents and ancestors in the coalescing loop :
  ancestorsDemes <- descendentDemes <- tipDemes 
  # A list of demes with all the genes remaining in each deme. Initialized and optimized.
  nodes_remaining_by_deme = list()
  deme <- as.array(seq(from=1,to=length(K),by=1))
  nodes_remaining_by_deme <- lapply(X=deme, FUN=remainingNodes, descendentDemes)
  
  # Number of coalescence events :
  single_coalescence_events=0 
  single_and_multiple_coalescence_events=0
  
  ### Simulating the coalescent process :
  while (length(unlist(nodes_remaining_by_deme))>1) 
  {
    
    ## Migration
    # we localize the ancestors in the landscape by sampling in the backward transition matrix. Optimized.
    names_node <- names(ancestorsDemes)
    ancestorsDemes <- apply(X = as.array(1:length(ancestorsDemes)), 
                            MARGIN=1, 
                            FUN = backwardParentsLocalizationSampling, 
                            rasterStack = rasterStack, 
                            transitionBackward = transitionBackward, 
                            DemeIdOfNodes = descendentDemes, 
                            K = K) 
    
    names(ancestorsDemes) <- names_node
    # once we know the ancestor deme numbers, we calculate the forward dispersion probability of the event
    time=time+1; 
    prob_forward[time] = sum(log(transitionForward[ancestorsDemes,descendentDemes]))
    number_of_nodes_over_generations = number_of_nodes_over_generations + length(descendentDemes)
    
    ## Coalescence
    
    # we now perform coalescence within each deme of the landscape for the ancestors
    for (deme in 1:length(K))#deme=1;deme=2;deme=3;deme=4;deme=5;deme=26;deme=10
    {
      # add a local variable, easier to manipulate than the reference to a list...
      nodes_remaining_in_the_deme = nodes_remaining_by_deme[[deme]] <- names(which(ancestorsDemes==deme))
      
      # we obtain the identities in the geneticData table (line) of the nodes remaining in the deme
      if (length(nodes_remaining_in_the_deme)>1)
      {
        # Create a function for ancestor attribution within a deme :
        ancestorDescendentmatrix <- parentalityAttributationWithinADeme(nodesRemainingInDeme = nodes_remaining_in_the_deme, N = round(K), deme=deme)
        
        # Columns of ancestorDescendentmatrix with more than one TRUE allow to identify coalescing individuals :
        if (any(colSums(ancestorDescendentmatrix)>1) )
        {
          #  Loop over all the ancestors in which coalescence event occur
          for (multiple in which(colSums(ancestorDescendentmatrix)>1)) # multiple<-which(colSums(ancestorDescendentmatrix)>1)[1]
          {
            # Record the coalescence event : 
            single_coalescence_events = single_coalescence_events +1
            
            # which(ancestorDescendentmatrix[,multiple]) identifies which node in the column coalesce
            nodes_that_coalesce = names(which(ancestorDescendentmatrix[,multiple]))
            
            # attibutes new node number to the ancestor
            new_node <- max(nodes)+1
            # removes the nodes that coalesced from the node vector
            nodes=nodes[!(names(nodes)%in%nodes_that_coalesce)]
            # adds them to the nodes vector
            nodes=append(nodes,new_node)
            names(nodes)[length(nodes)]=new_node
            
            # updating of vector ancestorsDemes (adding the deme number of the new node and removing the nodes that disapeared)
            ancestorsDemes <- append(ancestorsDemes[!(names(ancestorsDemes)%in%nodes_that_coalesce)],deme)
            names(ancestorsDemes)[length(ancestorsDemes)]<-new_node
            # adds the event to the list coalescent: time, which node coalesced, and the number of the new node
            coalescent[[single_coalescence_events]] <- list(time=time,coalescing=as.numeric(nodes_that_coalesce),new_node=new_node)
            # updating the nodes vector for the deme
            nodes_remaining_in_the_deme = nodes_remaining_by_deme[[deme]] <- append(nodes_remaining_in_the_deme[!nodes_remaining_in_the_deme %in% nodes_that_coalesce],new_node)
            # updates the number of coalescent events 
            single_and_multiple_coalescence_events = single_and_multiple_coalescence_events + length(nodes_that_coalesce) - 1
            
          } #  end of loop over all the ancestors in which coalescence event occur
        } # end of the if condition "there are coalescing events"
      } # end of the condition "there are more than 1 individual in the deme
    } # end of the loop across the demes
    
    descendentDemes = ancestorsDemes
  } # end of the backward generation while coalescence loop
  list(coalescent=coalescent,forward_log_prob=sum(prob_forward)/number_of_nodes_over_generations)
  # forward_log_prob is the average per generation of the log probability of the forward movements of the genes in the landscape
}


parentalityAttributationWithinADeme <- function(nodesRemainingInDeme, N, deme)
{
  # Create a function for parentality attribution within a deme used in simul_coalescent
  #
  # Args:
  #   nodesRemainingInDeme:
  #   N: a vector giving the population size in each deme
  #   deme: numeric giving the deme we want to analyse
  #
  # Returns:
  #   A matrix of parentality attribution. Two nodes coalesce if they have TRUE for the same parent (parents are in columns)
  
  nbGenesRemaining=length(nodesRemainingInDeme)
  
  # Attribute parents (among K possible parents) to each node present in the deme
  smp = sample(N[deme],length(nodesRemainingInDeme),replace=TRUE)
  
  # A logical matrix in which lines represent the nodes in the deme and column represent their parent :
  # (actually, this line is a simple test to transform the parentality info under a TRUE/FALSE form)
  parentDescendentmatrix <- matrix(smp, nrow=nbGenesRemaining, ncol=N[deme]) == matrix(1:N[deme], nrow=nbGenesRemaining, ncol=N[deme], byrow=TRUE)
  rownames(parentDescendentmatrix) <- nodesRemainingInDeme
  
  return(parentDescendentmatrix)
}

backwardParentsLocalizationSampling <- function(node, rasterStack, transitionBackward, DemeIdOfNodes, K)
{
  # Localizes the parents in the landscape by sampling in the backward transition matrix.
  #
  # Args:
  #   node: the node from which to set parents localization
  #   rasterStack: the rasterStack used to define landscape structure (demes)
  #   transitionBackward: the transition matrice giving the probability the migrate between demes of the landscape
  #   DemeIdOfNodes: the vector giving nodes localization (giving the number of the deme where each node sits)
  #   K: the vector giving the values of carrying capacities
  #
  # Returns:
  #   The number of the deme where the parents are likely to be found
  return(sample(length(K),size=1,prob=c(transitionBackward[DemeIdOfNodes[node],])))
}

remainingNodes <- function(deme, demeIdOfNodes)
{
  # Finds all the nodes remaining a specified deme
  #
  # Args:
  #   deme: a numerci giving the deme in which we want to find the remaining nodes
  #   demeIdOfNodes: gives the nodes localization
  return(which(demeIdOfNodes==deme))
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


add_br_length_and_mutation <- function(coalescent,mutation_rate)
{
  # Adds br_length and mutation to coalescent list
  # arguments : 
  # - coalescent: the coalescent list (ordered from present to past)
  # - mutation_rate: the number of mutation events per generation
  # - initial_genetic_value
  # value : the coalescent list with branch lengt and number of mutation events
  # 
  # Example
  # trB = matrix(c(1/4,1/2,1/4,0,1/3,1/3,1/6,1/6,1/2,1/4,1/8,1/8,1/5,2/5,2/5,0),nrow=4,ncol=4,byrow=TRUE)
  # trF = matrix(c(1/4,1/2,1/2,1/8,1/3,1/3,1/6,1/6,1/2,1/4,1/8,1/8,0,0,2/5,0),nrow=4,ncol=4,byrow=TRUE)
  # K=c(4,3,1,5)
  # tipsDemes = c(1,4,2,2,1,1,2,3);names(tipsDemes)=1:8
  # Coalescent = simul_coalescent_only(tipDemes=tipsDemes,transitionForward=trF,transitionBackward=trB,K=K)
  # Coalescent_genetics <- add_br_length_and_mutation(Coalescent,mutation_rate=.1)
  #
  #
  tips = NULL
  internals = NULL
  nodes = NULL
  times = NULL
  for (i in 1:length(coalescent[[1]]))#i=1;i=2
  {
    nodes = append(nodes,c(coalescent[[1]][[i]]$coalescing,coalescent[[1]][[i]]$new_node))
    internals = append(internals,coalescent[[1]][[i]]$new_node)
    times = append(times,coalescent[[1]][[i]]$time)
  }
  nodes = as.numeric(levels(as.factor(c(nodes,internals))));nodes = nodes[order(nodes)]
  tips = nodes[!((nodes)%in%(internals))]
  # getting the branch length of each coalescing node
  for (i in 1:length(coalescent[[1]]))#i=1
  {
    for (coalescing in coalescent[[1]][[i]]$coalescing)# coalescing = coalescent[[1]][[i]]$coalescing[1]
    {
      if (coalescing %in% tips) {coalescent[[1]][[i]]$br_length <- append(coalescent[[1]][[i]]$br_length,coalescent[[1]][[i]]$time)
      } else {
        coalescent[[1]][[i]]$br_length <- append(coalescent[[1]][[i]]$br_length,coalescent[[1]][[i]]$time-times[which(internals==coalescing)]) 
      } 
      coalescent[[1]][[i]]$mutations <- rpois(rep(1,length(coalescent[[1]][[i]]$br_length)),coalescent[[1]][[i]]$br_length*mutation_rate)
    }
  }
coalescent
}

coalist_2_coaltable <- function(coalist)
{
  # Conversion from coalist to coalecent table
  # Argument:
  # coalist : a list of coalescent events from present to past
  # with sublist "coalescing"= descendents nodes numbers, "new_node" = ancestor node number, br_length, mutations
  # Values:
  # coaltable : a table with the same information on coalescent events from present to past as lines
  #
  # example:
  # coalescent=list(list(time=1,coalescing=1:2,new_node=4,br_length=c(1,1),mutations = c(1,0)),list(time=4,coalescing=c(3:4),new_node=5,br_length=c(4,3),mutations = c(2,0)))
  # coalist_2_coaltable(coalescent)
  coaldf <- data.frame(Reduce(rbind,coalist))
  # note: warnings due to repeated lines names here
  coaltable <- coaldf[rep(1:dim(coaldf)[1],unlist(lapply(coaldf$coalescing, length))),]
  coaltable[,c("coalescing","br_length","mutations")] <- unlist(coaldf[,c("coalescing","br_length","mutations")])
  coaltable
}
