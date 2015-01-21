stepWiseMutationModel <- function(coaltable,initial_genetic_value,stepvalue)
{
  # Compute a change in the genetic value along a branch a a coalescent 
  #
  # Args: 
  #   coaltable: 
  #   initial_genetic_value:
  #   stepvalue:
  #
  # Returns:
  #   The coaltable
  coaltable$genetic_value=NA
  # we calculate the oritattion of the mutations in the different branches using binomial rules
  coaltable$resultant = 2*(rbinom(dim(coaltable)[1],coaltable[,"mutations"],.5)-coaltable[,"mutations"]/2)
  coaltable[dim(coaltable)[1]+1,] <- c(NA,max(unlist(coaltable$new_node)),NA,NA,NA,initial_genetic_value,NA)
  for(branch in rev(rownames(coaltable)[-dim(coaltable)[1]]))
  {
    coaltable[branch,"genetic_value"] <- coaltable[branch,"resultant"]*stepvalue + coaltable[which(coaltable$coalescing==coaltable[branch,"new_node"]),"genetic_value"]
  }
  coaltable
}

twoPhasesModel <- function(coaltable,initial_genetic_value,stepvalue,mut_param=c(p=.5,sigma2=4))
{
  p_loi_geometrique = ((1+4*mut_param["sigma2"])^.5-1)/(2*mut_param["sigma2"])
  coaltable$genetic_value=NA
  # we calculate the type of mutations in the different branches using binomial rules
  # either stepwise, or geometric positive, or geopetric negative
  coaltable[,c("n_stepw","n_geom_pos","n_geom_neg")] <- t(rmultinom(dim(coaltable)[1],coaltable[,"mutations"],c(mut_param["p"],(1-mut_param["p"])/2,(1-mut_param["p"])/2)))
  coaltable$resultant_stepw <- 2*(rbinom(dim(coaltable)[1],coaltable[,"n_stepw"],.5)-coaltable[,"n_stepw"]/2)
  coaltable$resultant_geom <- (rnbinom(dim(coaltable)[1],size=coaltable$n_geom_pos,p_loi_geometrique)
                               - rnbinom(dim(coaltable)[1],size=coaltable$n_geom_neg,p_loi_geometrique))
  coaltable$resultant_geom[is.na(coaltable$resultant_geom)]=0
  coaltable[dim(coaltable)[1]+1,] <- c(NA,max(unlist(coaltable$new_node)),NA,NA,NA,initial_genetic_value,NA,NA,NA,NA,NA)
  for(branch in rev(rownames(coaltable)[-dim(coaltable)[1]]))
  {
    coaltable[branch,"genetic_value"] <- (coaltable[branch,"resultant_stepw"]+coaltable[branch,"resultant_geom"])*stepvalue + coaltable[which(coaltable$coalescing==coaltable[branch,"new_node"]),"genetic_value"]
  }
  coaltable
}


bigeometricModel <- function(coaltable,initial_genetic_value,stepvalue,mut_param=c(sigma2=4))
{
  p_loi_geometrique = ((1+4*mut_param["sigma2"])^.5-1)/(2*mut_param["sigma2"])
  coaltable$genetic_value=NA
  # we calculate the type of mutations in the different branches using binomial rules
  # either stepwise, or geometric positive, or geopetric negative
  coaltable[,"n_geom_pos"] <- rbinom(dim(coaltable)[1],coaltable[,"mutations"],c(mut_param["p"],(1-mut_param["p"])/2,(1-mut_param["p"])/2))
  coaltable$resultant_geom <- (rnbinom(dim(coaltable)[1],size=coaltable$n_geom_pos,p_loi_geometrique)
                               - rnbinom(dim(coaltable)[1],size=coaltable$mutations-coaltable$n_geom_pos,p_loi_geometrique))
  coaltable$resultant_geom[is.na(coaltable$resultant_geom)]=0
  coaltable[dim(coaltable)[1]+1,] <- c(NA,max(unlist(coaltable$new_node)),NA,NA,NA,initial_genetic_value,NA,NA,NA,NA)
  for(branch in rev(rownames(coaltable)[-dim(coaltable)[1]]))
  {
    coaltable[branch,"genetic_value"] <- coaltable[branch,"resultant_geom"]*stepvalue + coaltable[which(coaltable$coalescing==coaltable[branch,"new_node"]),"genetic_value"]
  }
  coaltable  
}