stepWiseMutationModel <- function(mutations)
{
  # Compute a resultant according to step wise mutation model
  #
  # Args: 
  #   mutations: number of mutations
  #
  # Returns:
  #   The resultant in microsatellite repetition number using binomial rules
  
  res <- 2*rbinom(n = length(mutations), size = mutations, prob = 0.5) - mutations
  return(res)
  
  # Example:
  # stepWiseMutationModel(mutations=5)
  # stepWiseMutationModel(mutations=c(2,3))
}

twoPhasesModel <- function(mutations=3,p=.5,sigma2=4)
{
  # Simulates a change in the genetic value depending on number of mutations
  # assuming two phases mutation model
  #
  # Args: 
  #   mutations= number of mutations
  #
  # Value
  #   The resultant in microsatellite repetition number using binomial rules
  # 
  # Example:
  # twoPhasesModel(5)
  #
  # Estimate parameter of geometric mutation size
  p_loi_geometrique = ((1+4*sigma2)^.5-1)/(2*sigma2)
  # Sample number of stepwise mutations among mutations from binomial law
  NStepw <- rbinom(length(mutations),mutations,prob=c(p,1-p))
  # Estimate pascal triangle resultant of stepwise mutations 
  resultantStepw <- 2*rbinom(n = length(mutations), size = NStepw, prob = 0.5) - NStepw
  # Estimate number of geometric mutations having positive effect on microsatallite size
  NGeomPos <- rbinom(length(mutations),mutations-NStepw,.5)
  # Estimate resultant geometric mutation changes for positive and negative effects
  resultantGeom <- rnbinom(length(mutations)*2,size=c(NGeomPos,mutations-NStepw-NGeomPos),p_loi_geometrique)
  # note: warnings due to size=0 in rnbinom parameters
  resultantGeom[is.na(resultantGeom)] <- 0
  # Overall resultant is binomial plus geometric positive minus geometric negative
  return(resultantStepw+c(diff(t(matrix(resultantGeom,nrow=length(mutations),ncol=2)))))
}


bigeometricModel <- function(mutations=3,sigma2=4)
{
  # Simulates a change in the genetic value depending on number of mutations
  # assuming geopetric model
  #
  # Args: 
  #   mutations= number of mutations
  #
  # Value
  #   The resultant in microsatellite repetition number using binomial rules
  # 
  # Example:
  # bigeometricModel(5)
  #
  # Estimate parameter of geometric mutation size
  p_loi_geometrique = ((1+4*sigma2)^.5-1)/(2*sigma2)
  # Estimate number of geometric mutations having positive effect on microsatallite size
  NGeomPos <- rbinom(length(mutations),mutations,.5)
  # Estimate resultant geometric mutation changes for positive and negative effects
  resultantGeom <- rnbinom(length(mutations)*2,size=c(NGeomPos,mutations-NGeomPos),p_loi_geometrique)
  # note: warnings due to size=0 in rnbinom parameters
  resultantGeom[is.na(resultantGeom)] <- 0
  # Overall resultant is binomial plus geometric positive minus geometric negative
  return(diff(t(matrix(resultantGeom,nrow=length(mutations),ncol=2))))
}


resultantFunction <- function(nbrMutations, stepValue, mutationModel, args){
  # Compute a resultant giving a vector of number of mutation
  #
  # Args :
  #   nbrMutations : a vector giving the number of mutations
  #   stepValue : the stepValue of the locus
  #   mutationModel : the mutation model to apply
  #   args : a list containing the mutation models parameters values
  #
  # Returns:
  #   a vector of resultant
  res <- do.call(what = mutationModel, args = c(list(nbrMutations), args))
  res <- res * stepValue
  return(res)
}

addGeneticValueToCoaltable <- function(coalTable,initialGenetValue,stepValue)
{
  # Compute genetic value of nodes in a coalescent, starting from the ancestor
  #
  # Args:
  #   coalTable: a table giving the coalescent, returned by coalist_2_coaltable function
  #   initialGenetValue: the genetic value of the ancestor
  #   stepValue : the step value of the mutation model
  #
  # Returns:
  #   the appended coalTable, with genetic value
  
  # Add a row for the ancestor
  coalTable[dim(coalTable)[1]+1,"genetic_value"] <- initialGenetValue
  # Precise the children nodes of the ancestor to connect the genealogy to the ancestral genetic value
  coalTable[dim(coalTable)[1],"coalescing"] <- max(unlist(coalTable[,"new_node"]))
  
  # Compute genetic values from ancestor to children
  for(branch in rev(rownames(coalTable)[-dim(coalTable)[1]])) # branch <- rev(rownames(coalTable)[-dim(coalTable)[1]])[1]
  {
    # for the child, add parent genetic value and resultant to get child genetic value
    coalTable[branch,"genetic_value"] <- coalTable[branch,"Resultant"] + coalTable[which(coalTable$coalescing==coalTable[branch,"new_node"]),"genetic_value"]
  }
  
  return(coalTable)
}
