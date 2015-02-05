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
  dist(geneticData)
}


PCA_rotation <- function(geneticData,Distance)
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
  GenetDist = do.call(Distance,geneticData)
  rotation = prcomp(GenetDist)$rotation
  return(rotation)
  # note we may have to select majors PCi containing information
}

