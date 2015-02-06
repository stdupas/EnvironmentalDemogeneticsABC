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
  p_loi_geometrique = ((1+4*sigma2)^.5-1)/(2*sigma2)
  n_stepwN_geomN_geom_neg <- t(rmultinom(length(mutations),mutations,c(p,(1-p)/2,(1-p)/2)))
  resultant_geom <- rnbinom(length(mutations)*2,size=c(N_geom_pos,mutations-N_geom_pos),c(p_loi_geometrique,p_loi_geometrique))  
  # note: warnings due to size=0 in rnbinom parameters
  resultant_geom[is.na(resultant_geom)] <- 0
  resultant_stepw+c(diff(t(matrix(resultant_geom,nrow=length(mutations),ncol=2))))
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
  p_loi_geometrique = ((1+4*sigma2)^.5-1)/(2*sigma2)
  N_geom_pos <- rbinom(length(mutations),mutations,.5)
  resultant_geom <- rnbinom(length(mutations)*2,size=c(N_geom_pos,mutations-N_geom_pos),c(p_loi_geometrique,p_loi_geometrique))  
  # note: warnings due to size=0 in rnbinom parameters
  resultant_geom[is.na(resultant_geom)] <- 0
  c(diff(t(matrix(resultant_geom,nrow=length(mutations),ncol=2))))
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

addGeneticValueToCoaltable <- function(coaltable,initialGenetValue,stepvalue)
{
  coaltable[dim(coaltable)[1]+1,"genetic_value"]=initialGenetValue
  coaltable[dim(coaltable)[1],"coalescing"]=max(unlist(coaltable[,"new_node"]))
  for(branch in rev(rownames(coaltable)[-dim(coaltable)[1]])) # branch=rev(rownames(coaltable)[-dim(coaltable)[1]])[1]
  {
    coaltable[branch,"genetic_value"] <- coaltable[branch,"Resultant"] + coaltable[which(coaltable$coalescing==coaltable[branch,"new_node"]),"genetic_value"]
  }
coaltable
}
