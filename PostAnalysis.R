
pca4abc <- function(GeneticData, ParamList, distanceMethod, path){
  # From real and simulated data, computes via PCA necessary elements for further analysis
  #
  # Args:
  #   GeneticData : the genetic data set, individuals in rows and in columns x,y,Locus1... LocusN
  #   ParamList: the list of models parameters values get with the askListOfParameters function
  #   distanceMethod
  #
  # Return:
  #   a list of a vector of observed statistics, an array of simulated summary stats and a vector of parameters values
  
  
  # number of loci under study:
  obsGenetics <- GeneticData[!(colnames(GeneticData)%in%c("x","y","Cell_numbers"))]
  
  # number of individuals under study:
  nbrInd <- nrow(obsGenetics)
  
  ### Caracterizing rotation to apply to genetic simulations to get summary statistics
  rotation = PCA_rotation(geneticData = obsGenetics, DistanceMethod = distanceMethod)
  
  summaryStatObsMat = as.matrix(do.call(what = distanceMethod, args = list(obsGenetics)))%*%rotation
  
  colN <- paste(rep(colnames(summaryStatObsMat), each = length(rownames(summaryStatObsMat))),
                rownames(summaryStatObsMat),sep=".")
  
  summaryStatObsVect = c(summaryStatObsMat)
  
  names(summaryStatObsVect)<- colN
  
  ### Computing summary statistics of simulated data
  allFiles <- grep(pattern = "^Genetics_\\d*.txt$", x=list.files(path), value = TRUE)
  
  stats <- apply(X = as.array(allFiles), MARGIN = 1, FUN = computeSummaryStats, 
                 nbrInd = nbrInd, 
                 distanceMethod = distanceMethod, 
                 rotation = rotation,
                 path = path)

  stats <- stats[, order(stats[1,])]
  
  
  ### Return necessary elements for abc analysis :
  
  # First element for package abc
  statobs <- as.data.frame(t(summaryStatObsVect))
  
  # Second element for package abc
  sumstat <- t(stats)[,-1]
  colnames(sumstat) <- names(statobs)
  
  # Third element for package abc :
  temp <- as.data.frame(ParamList)
  param <- temp[, grep(pattern = "Values", x = names(temp))]
  
  return(list(statobs = statobs, sumstat = sumstat, param = param))
}

