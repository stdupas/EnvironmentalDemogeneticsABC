SharedAlleleDistance <- function(geneticData)
{
  # shared allele genetic distance function (Chakraborty and Jin, 1993)
  # calculates proportion of shared alleles between haploid individuals among locus
  # argument : geneticData (with columns as loccus named "LocusX" and lines as haplotypes)
  geneticDataArray <- array(unlist(c(geneticData)), dim=c(dim(geneticData)[1],dim(geneticData)[2],dim(geneticData[1])))[,,,1]
  dist(rowMeans(aperm(geneticDataArray==aperm(geneticDataArray,c(3,2,1)),c(1,3,2)),dims=2))
}

DeltaMuDistance <- function(geneticData)
{
  dist(geneticData)^2 / dim(geneticData)[2]
}


PCA_rotation <- function(geneticData,DistanceMethod)
{
  # Calculates the rotation to apply to genetic data for summary stats 
  # argument:
  # geneticDataObs: observed genetic data (columns x, y, Cell_numbers, Locus1 ... Locusn)
  # value: 
  # rotation corresponding to the PCA axis
  #
  # Example:
  # GeneticData = data.frame(x=c(1.5,1.5,3.5,3.5),y=c(.5,.5,1.5,1.5),Locus1=c(200,204,206,206),Locus2=c(156,154,166,164))
  # PCA_rotation(GeneticData,DeltaMuDistance)
  # PCA_rotation(GeneticData,SharedAlleleDistance)
  genetDist = do.call(what = DistanceMethod, args = list(geneticData))
  PrincipalComponentResults = prcomp(genetDist)
  lastsignificant = 1+max(which(summary(PrincipalComponentResults)$importance["Cumulative Proportion",]<.95))
  rotation = PrincipalComponentResults$rotation[,1:lastsignificant]
  return(rotation)
  # note we may have to select majors PCi containing information
}


appendGenetRefTable <- function(simulation=1,rotation,locusNames,DistanceMethod)
{
  # adds summary statistics for individual i
  # 
  # arguments
  # 
  # 
  #
  geneticData <- read.table(paste("SimulResults/","Genetics_",simulation,".txt",sep=""))[,locusNames]
  GenetDist = do.call(what = DistanceMethod, args = list(geneticData))
  
}

computeSummaryStats <- function(file, nbrInd, distanceMethod, rotation){
  # Compute the summary statistics of a genetic file (with the forward log likelyhood)
  #
  # Args : 
  #   file: the file to be read, containing individuals by row and loci by columns, and one line of genetic value
  #   nbrInd: the number of individuals under study
  #   distanceMethod : the method to be applied to compute genetic distances
  #   rotation: the pca rotation to be applied to genetic distances
  #
  # Returns: 
  #   a vector containing : the ID of the simulation file, the forwardLogLikelyhood, the projected distances
  
  # Read file
  genetics <- read.table(file = file, nrows = nbrInd)
  
  # Compute distance
  genetDist = do.call(what = distanceMethod, args = list(genetics))
  rotatedGenetDist = as.matrix(genetDist)%*%rotation
  
  # Extract the exact number of the simulation
  sim <- as.numeric(gsub(pattern = "[^0-9]", replacement = "", x = file))
  return(c(sim, rotatedGenetDist))
}