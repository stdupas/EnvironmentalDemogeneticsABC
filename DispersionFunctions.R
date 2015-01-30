################### Formated Dispersion Functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fatTail1 <- function(x, alpha, beta){
  # Computes a value for kernel dispersion using a fat tail model, according to the method in ?
  #
  # Args: 
  #   x: the distance between two points
  #   alpha : the first parameter of the dispersion law
  #   beta: the second parameter of the dispersion law
  #
  # Returns:
  #   The value of dispersion kernel for x
  return(1/(1+1/alpha *x^beta))
}

fatTail2 <- function(x, sigma, gamma){
  # Computes a value for kernel dispersion using a fat tail model, according to the method in ?
  #
  # Args: 
  #   x: the distance between two points
  #   sigma : the first parameter of the dispersion law
  #   gamma: the second parameter of the dispersion law
  #
  # Returns:
  #   The value of dispersion kernel for x
  return(x^gamma*exp(-2*x/(sigma^0.5)))
}

gaussian <- function(x, sigma){
  # Computes a value for kernel dispersion using a gaussian model, i.e. a simple normal density distribution (sigma, mean=0)
  #
  # Args: 
  #   x: the distance between two points
  #   sigma: the value of the standard deviation
  #
  # Returns:
  #   The value of dispersion kernel for x
  return(dnorm(x, mean = 0, sd = sigma, log = FALSE))
}

exponential <- function(x, sigma){
  # Computes a value for kernel dispersion using an exponential model, at rate=1/sigma
  #
  # Args: 
  #   x: the distance between two points
  #   sigma: the value of the standard deviation
  #
  # Returns:
  #   The value of dispersion kernel for x
  
  return(dexp(x, rate = 1/sigma, log = FALSE))
}

contiguous <- function(x, sigma, threshold){
  # Computes a value for kernel dispersion using a contiguous model (near dispersal)
  #
  # Args: 
  #   x: the distance between two points
  #   sigma: ?
  #   threshold : a trashold value
  #
  # Returns:
  #   The value of dispersion kernel for x
  return((x==0)*(1-sigma)+((x>0)-(x>1.4*threshold))*(sigma/2))
}

island <- function(x, sigma){
  # Computes a value for kernel dispersion using an island model (probability 1-m to stay, else homogen dispersion)
  #
  # Args: 
  #   x: the distance between two points
  #   sigma: ?
  #
  # Returns:
  #   The value of dispersion kernel for x
  return((x==0)*(1-sigma)+(x>0)*(sigma))
}

gaussianMixedIsland <- function(x, sigma, gamma, normalizingFactor){
  # Computes a value for kernel dispersion using a mixed model (gaussian/island)
  #
  # Args: 
  #   x: the distance between two points
  #   sigma: the standard deviation of the gaussian (normal) function
  #   gamma: the parameter of the island model
  #   normalizaingFactor: a factor in practice equal to the product of the migration matrix dimensions)
  # Returns:
  #   The value of dispersion kernel for x
  (gamma+(1-gamma)*dnorm(x, mean = 0, sd = sigma, log = FALSE))/prod(dim(migration))
}

################# End of formated Dispersion Functions <<<<<<<<<<<<<<<<<<<<


################# Applying Dispersion Functions to Objects >>>>>>>>>>>>>>>>>>
dispersionFunctionForValue <- function(dispersionFunction, x, args){
  # Compute a dispersion kernel function over a single value
  #
  # Args:
  #   dispersionFunction: the name of the dispersion kernel function which is called
  #   x: the distance between two points
  #   args : a list of the arguments of the dispersion function
  # 
  # Returns:
  #   The value of dispersion kernel for x
  args <- c(list(x), args)
  return(do.call(dispersionFunction, args))
  
  # Ex : 
  # dispersionFunctionForValue(fatTail1, x=4, args=list(alpha=0.5, beta=0.2))
}

dispersionFunctionForArray <- function(dispersionFunction, array, args){
  # Apply a dispersion kernel function over an array
  #
  # Args:
  #   dispersionFunction: the name of the dispersion kernel function which is called
  #   array: the array of distances used to compute dispersion function.
  #   args : a list of the arguments of the dispersion kernel function
  # 
  # Returns:
  #   An array corresponding to the dispersion kernel values
  apply(X=array, MARGIN=1, FUN=dispersionFunctionForValue, args=args, dispersionFunction=dispersionFunction)
  
  # Ex:
  # dispersionFunctionForArray(dispersionFunction=fatTail2,
  #                            array=array(data= 1:10, dim =10),
  #                            args=list(sigma = 0.2, gamma=0.3 ))
}

dispersionFunctionForRasterLayer <- function(dispersionFunction, rasterLayer, args){
  # Apply a dispersion kernel function over distances between the cells of a rasterLayer
  #
  # Args:
  #   dispersionFunction: the dispersion kernel function to apply
  #   rasterLayer: a rasterLayer with coordinates set
  #   args: a list of the dispersion function arguments
  #
  # Returns:
  #   A matrix of dispersion kernel values, similar in shape with the distance matrix between cells.
  
  # Extract coordinates from rasterLayer
  coords = xyFromCell(object=rasterLayer, cell=1:length(values(rasterLayer[[1]])), spatial=FALSE)
  # Compute distance matrix
  distanceMatrix = as.matrix(dist(coords))
  # Apply dispersionFunction
  dispersion <- apply(X=distanceMatrix, 
                      MARGIN=c(1,2), 
                      FUN=dispersionFunctionForValue, 
                      dispersionFunction=dispersionFunction, 
                      args=args)
  return(dispersion)
}

migrationRateMatrix <- function(dispersion){
  # Normalizes a matrix of dispersion kernel between cells to get a migration rate matrix between cells.
  #
  # Args:
  #   dispersion: a matrix representing the values of a specified kernel (function of distances between cells)
  #
  # Returns:
  #   A migration rate matrix (note that rowSums and colSums are not 1: cause of bordure effect, individuals go "out of the world")
  return(dispersion/max(c(colSums(dispersion),rowSums(dispersion))))
}