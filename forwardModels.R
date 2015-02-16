forwardMigrationRateMatrixFromKernel <- function(dispersionKernel){
  # Normalizes a matrix of dispersion kernel between cells so that individuals bounce on the frontier.
  #
  # Args:
  #   dispersion: a matrix representing the values of a specified kernel (function of distances between cells)
  #
  # Returns:
  #   A migration rate matrix (note that rowSums and colSums are not 1: cause of bordure effect, individuals go "out of the world")
  # ex: 
  # kernelMatrix=c(1.000000e+00,0.804259048,0.004416269,8.156667e-05,
  #                 8.042590e-01,1.000000000,0.804259048,4.416269e-03,
  #                 4.416269e-03,0.804259048,1.000000000,8.042590e-01,
  #                 8.156667e-05,0.004416269,0.804259048,1.000000e+00),nrow=4,ncol=4)
  # fmatrix <- forwardMigrationRateMatrixFromKernel(kernelMatrix)
  # rowSums(fmatrix)
  return(dispersionKernel/matrix(rowSums(dispersionKernel),nrow=nrow(dispersionKernel),ncol=ncol(dispersionKernel),byrow=FALSE))
}

reproduction <- function(tipDemes, r)
{
  # Reproduces individuals according to growth rate of their deme
  # Arguments:
  # - tipeDemes : a data frame describing the individuals to reproduce
  # with columns "individualNb" (individual number) and "demeNb" (deme number)
  # - r : a vector of growth rate of each deme
  # Value:
  # A data.frame describing the offsprings, 
  # - column "individualNb" gives the number of the parent, 
  # - column "demeNb" gives the deme 
  # rownames have the format "parent_number.offspring_number"
  # 
  # Example
  # r<-2^(0:3)
  # tipDemes=data.frame(individualNb=1:10,demeNb=sample(1:4,10,replace=TRUE))
  # reproduction(tipDemes,r)
  rownames(tipDemes)<-tipDemes$individualNb
  new_individuals <- rep(tipDemes$individualNb,rpois(nrow(tipDemes),r[tipDemes$demeNb ]))
  newTipDemes <- tipDemes[new_individuals,]
  return(newTipDemes)
}

absoluteForwardTransition <- function(r,K,migrationMatrix)
{
  # Calculte absolute forward transition matrix from r, K and forward migration matrix
  # 'absolute' means that rowSums is not 1, it determines positions of the offspring
  # and offspring population size can be different from that of the parent
  #
  # Arguments:
  # - tipeDemes : a data frame describing the individuals to reproduce
  # with columns "individualNb" (individual number) and "demeNb" (deme number)
  # - r : a vector of growth rate of each deme
  # - K : a vector of carrying capacity of each deme
  # Value:
  # A transiton matrix describing parent offspring reproduciton and movements
  # 
  # Example
  # r<-2^(0:3)
  # K<-c(10,5,20,2)
  # migrationMatrix <- matrix(c(.9,.1,0,0,.05,.9,0.05,0,0,.05,.9,0.05,0,0,.1,.9),nrow=4,ncol=4,byrow=TRUE)
  # absoluteTransitionForward(r,K,migrationMatrix)
  #
  #
  absTrF <- matrix(0,nrow=length(r),ncol=length(r))
  sum_k.m_kj <- rep(0,length(r))
  for (k in 1:length(r))
    for (j in 1:length(r))
      sum_k.m_kj[j] <- sum_k.m_kj[j] + r[k]*migrationMatrix[k,j]
  for (i in 1:length(r))
    for (j in 1:length(r))
    {
      absTrF[i,j] <- r[i]*migrationMatrix[i,j]*min(1,K[j]/(sum_k.m_kj[j]))
    }
return(absTrF)
}

reprMigrCompet <- function(tipDemes,absoluteForwardTransitionMatrix,generationTime,generationTimeSD)
{
  # Randomly reproduces and migrates individuals according to growth rate of their deme
  # migration rates to descendent demes and carrying capcacity of descendent dames
  # Arguments:
  # - tipeDemes : a data frame describing the individuals to reproduce
  # with columns "individualNb" (individual number), "demeNb" (deme number) and "birthDate"
  # - absoluteForwardTransitionMatrix : absolute transition matrix with expected offspring size in column for each deme in line
  # - generationTime : a vector of generation time for each deme
  # - generationTimeRelativeSD : a standard error relative to the mean for generation time
  # Value:
  # A data.frame describing the offsprings, 
  # - column "individualNb" gives the number of the parent, 
  # - column "demeNb" gives the deme 
  # - column "birthdate"
  # rownames have the format "parent_number.offspring_number"
  # 
  # Example
  # r<-2^(0:3)
  # K<-c(10,5,20,2)
  # tipDemes=data.frame(individualNb=1:10,demeNb=sample(1:4,10,replace=TRUE),birthDate=as.Date("2013/01/23"))
  # generationTime=sample(20:25,4,replace=TRUE)
  # generationTimeRelativeSD=.1
  # migrationMatrix <- matrix(c(.9,.1,0,0,.05,.9,0.05,0,0,.05,.9,0.05,0,0,.1,.9),nrow=4,ncol=4,byrow=TRUE) 
  # reprMigrCompet(tipDemes,
  #                absoluteForwardTransitionMatrix=absoluteForwardTransition(r,K,migrationMatrix),
  #                generationTime,generationTimeRelativeSD)
  #

  rownames(tipDemes)<-tipDemes$individualNb
  # caculating absolute transition matrix with expected offsrpog size in column for each deme of the tip in line
  trTip <- absoluteForwardTransitionMatrix[tipDemes$demeNb,]
  # generating number of offspring of each tip (line) in each deme (column)
  individualNb <- matrix( rpois(prod(dim(trTip)),trTip),nrow=nrow(trTip),ncol=ncol(trTip))
  # 
  offsTipDemes <- tipDemes[rep(tipDemes$individualNb,rowSums(individualNb)),]
  offsTipDemes$birthDate <- offsTipDemes$birthDate + rnorm(nrow(offsTipDemes),generationTime[offsTipDemes$demeNb],generationTimeRelativeSD*generationTime[offsTipDemes$demeNb])
  offsTipDemes$demeNb <- rep(rep(1:length(r),nrow(individualNb)),as.vector(t(individualNb)))
  return(offsTipDemes)
}

generate_parameterSerie <- function(rstSeries,startingDate,stoppingDate,periodicity)
{
  
}

generateIndividualSerie <- function(release,r_serie,K_serie,migrationMatrix, startingDate, stoppingDate,generationTime,generationTimeSD)
{
  # function to generate a time serie of individuals birth date and deme
  # arguments:
  # r_serie : matrix of growh rate per deme in line and day in column
  # K_serie : matrix of carrying capacity per deme in line and day in column
  # release : a data frame describing the individuals to reproduce
  # startingDate, stoppingDate: dates to start and stop the simulation
  # stoppingDate: date to stop the time serie simulation
  # with columns "individualNb" (individual number), "demeNb" (deme number) and "birthDate"
  # value:
  # a data frame with all individuals characteristics in the time serie
  # example:
  #
  # r_serie = matrix(rpois(731*4,c(3,4,5,7)),
  #                  nrow = 4, ncol = 731, byrow=TRUE,
  #                  dimnames = list(1:4,as.character(as.Date(as.Date("2001/01/01"):as.Date("2003/01/01"),origin="1970-01-01"))))
  # K_serie = matrix(10,nrow = 4, ncol = 731, byrow=TRUE,
  #                  dimnames = list(1:4,as.character(as.Date(as.Date("2001/01/01"):as.Date("2003/01/01"),origin="1970-01-01"))))
  # release = data.frame(individualNb=1:10,demeNb=sample(1:4,10,replace=TRUE),
  #                     birthDate=as.Date(sample(as.Date(as.Date("2001/01/01"):as.Date("2001/02/01"),
  #                                       origin="1970-01-01"),10,replace=TRUE)))
  # startingDate=as.Date("2001/01/01")
  # stoppingDate=as.Date("2003/01/01")
  # generationTime=21,generationTimeSD=3
  individuals <- release
  for (Date in as.Date(startingDate:stoppingDate,origin="1970-01-01"))
  {
    absoluteForwardTransitionMatrix <- absoluteForwardTransition(r_serie[,as.character(Date)],K_serie[,as.character(Date)],migrationMatrix)
    individuals <- rbind(individuals,reprMigrCompet(individuals[which(individuals$birthDate==Date),],absoluteForwardTransitionMatrix,generationTime,generationTimeSD))
  }
individuals
}
