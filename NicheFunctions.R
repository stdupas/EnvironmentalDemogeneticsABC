################### Formated Niche Functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
binaryMultiplicative <- function(x, Y){
  # Norm reaction for binary environment variable used in a multiplicative response framework
  #
  # Args:
  #   x: binary value of the environmental variable : 0 or 1
  #   Y : value of the reaction norm for x = 1
  #
  # Returns: 
  #   The value of the reaction norm for
  res <- x*(Y-1) +1
  return(res)
}

conquadraticSkewed1 <- function(x, Xmin, Xmax, Xopt, Yopt)
{
  # Asymetric concave conquadratic function within an enveloppe, else returns 0.
  # 
  # Args: 
  #   x : numeric providing the values of variable to calculate reaction norm
  #   Xmin: start of the enveloppe
  #   Xmax: end of the enveloppe
  #   Xopt: the value that maximises the function
  #   Yopt: the maximum value of the function
  #
  # Returns:
  #   The value of reaction norm
  alpha <- -log(2)/log((Xopt-Xmin)/(Xmax-Xmin))
  Xprime<- ((x-Xmin)/(Xmax-Xmin))^alpha*(Xmax-Xmin)+Xmin
  y <- (Yopt-4*Yopt/((Xmax-Xmin)^2)*(Xprime-(Xmin+Xmax)/2)^2)*((x>=Xmin)&(x<=Xmax))
  y[x<Xmin] <-0
  y
}

conquadraticSkewed2 <- function(x, Xmin, Xmax, Xopt, Yopt)
{
  # Asymetric concave conquadratic function within an enveloppe, else returns 0.
  # 
  # Args: 
  #   x : numeric providing the values of variable to calculate reaction norm
  #   Xmin: start of the enveloppe
  #   Xmax: end of the enveloppe
  #   Xopt: the value that maximises the function
  #   Yopt: the maximum value of the function
  #
  # Returns:
  #   The value of reaction norm
  alpha <- -log(2)/log((Xopt-Xmin)/(Xmax-Xmin))
  Xprime<- ((x-Xmin)/(Xmax-Xmin))^alpha*(Xmax-Xmin)+Xmin
  y <- (Yopt-4*Yopt/((Xmax-Xmin)^2)*(Xprime-(Xmin+Xmax)/2)^2)*((x>=Xmin)&(x<=Xmax))
  y+y*(Yopt-y)
  #  y[X<Xmin] <-0
}


concaveQuadratic <- function(x, Xmin, Xmax, Yopt)
{
  # Concave quadratic function within an enveloppe, else returns 0.
  #
  # Args: 
  #   x: numeric providing the values of variable to calculate reaction norm
  #   Xmin: start of the enveloppe
  #   Xmax: end of the enveloppe
  #   Yopt: the maximum value of the function
  #
  # Returns:
  #   The value of the reaction norm
  (Yopt-(4*Yopt/(Xmax-Xmin)^2)*(x-(Xmin+Xmax)/2)^2)*((x>Xmin)&(x<Xmax))
}


concaveSquaredQuadratic <- function(x, Xmin, Xmax, Yopt)
{
  # Concave squared quadratic function within an eveloppe, else returns 0.
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   Xmin: start of the enveloppe
  #   Xmax: end of the enveloppe
  #   Yopt: the maximum value of the function
  #
  # Returns:
  #   The value of the reaction norm
  y = (Yopt-(4*Yopt/(Xmax-Xmin)^2)*(x-(Xmin+Xmax)/2)^2)*((x>Xmin)&(x<Xmax))
  y + y * (1-y)
}


AllOrNothing <- function(x, Xmin, Xmax, Yopt)
{
  # Returns Yopt if x is between Xmin and Xmax, else returns 0.
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   Xmin: start of the enveloppe
  #   Xmax: end of the enveloppe
  #   Yopt: the maximum value of the function
  #
  # Returns:
  #   The value of the reaction norm
  Yopt*((x>Xmin)&x<Xmax)
}


linearFourParameters <- function(x, Xmin, Xmax, Yxmin, Yxmax)
{
  # Computes a linear response within an enveloppe, else returns 0.
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   Xmin: start of the enveloppe
  #   Xmax: end of the enveloppe
  #   Yxmin: the value of the function at x=Xmin
  #   Yxmax: the value of the function at x=Xmax
  #
  # Returns:
  #   The value of the reaction norm
  a = (Yxmin - Yxmax) / (Xmin - Xmax)
  b = Yxmin - Xmin * a
  return((a*x+b)*((x>Xmin)&(x<Xmax)))
}

LogLinearFourParameters <- function(x, Xmin, Xmax, Yxmin, Yxmax)
{
  # Computes the logarithm of a linear response within an enveloppe, else returns 0.
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   Xmin: start of the enveloppe
  #   Xmax: end of the enveloppe
  #   Yxmin: the value of the function at x=Xmin
  #   Yxmax: the value of the function at x=Xmax
  #
  # Returns:
  #   The value of the reaction norm
  a = (Yxmin - Yxmax) / (Xmin - Xmax)
  b = Yxmin - Xmin * a
  return(log(a*x+b)*((x>Xmin)&(x<Xmax)))
}


linearTwoParameters <- function(x,X0,slope)
{
  # Computes a linear response within the enveloppe, else returns 0.
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   X0: value of the function at x=0
  #   slope: the value of the slope of the function
  #
  # Returns:
  #   The value of the reaction norm
  return(slope*x-slope*X0)
}

LogLinearTwoParameters <- function(x,X0,slope)
{
  # Computes a log linear response within the enveloppe, else returns 0.
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   X0: value of the function at x=0
  #   slope: the value of the slope of the function
  #
  # Returns:
  #   The value of the reaction norm
  return(log(slope*x-slope*X0))
}

linearPositiveTwoParameters <- function(x, X0, slope)
{
  # Compute a linear response . If response is negative, response is given value 0
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   X0: value of the function at x=0
  #   slope: the value of the slope of the function
  #
  # Returns:
  #   The value of the reaction norm
  return((slope*x-slope*X0)*(slope*x-slope*X0>0))
}

logLinearPositiveTwoParameters <- function(x, X0, slope)
{
  # Compute a log linear response . If response is negative, response is given value 0
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   X0: value of the function at x=0
  #   slope: the value of the slope of the function
  #
  # Returns:
  #   The value of the reaction norm
  return(log((slope*x-slope*X0)*(slope*x-slope*X0>1)+(slope*x-slope*X0<=1)))
}

constant <- function(x,Y)
{
  # Compute a constant response
  #
  # Args:
  #   x: numeric providing the values of variable to calculate reaction norm
  #   Y: value of the constant response
  #
  # Returns:
  #   The value of the reaction norm
  return(Y)
}

################# End of formated Niche Functions <<<<<<<<<<<<<<<<<<<<


################# Applying Niche Functions to Objects >>>>>>>>>>>>>>>>>>

nicheFunctionForValue <- function(nicheFunction, x, args){
  # Function to apply a niche function over a single value
  #
  # Args:
  #   nicheFunction: the name of the niche function which is called
  #   x: the value of the environemental variable
  #   args : a list of the arguments of the niche function
  # 
  # Returns:
  #   The value corresponding to the norm reaction for x
  args <- c(list(x), args)
  return(do.call(nicheFunction, args))
  
  # Ex : 
  # nicheFunctionForValue(conquadraticSkewed1, x=4, args=list(Xmin=0, Xmax=10, Xopt=5, Yopt=1))
}

nicheFunctionForArray <- function(nicheFunction, Array, args){ 
  # Function to apply a niche function over an Array
  #
  # Args:
  #   nicheFunction: the name of the niche function which is called
  #   array: the array of environmental values used to compute niche function.
  #   args : a list of the arguments of the niche function
  # 
  # Returns:
  #   An array corresponding to the norm reaction
  
  return(apply(X=Array, MARGIN=1, FUN=nicheFunctionForValue, nicheFunction=nicheFunction, args=args))
  # Ex:
  # nicheFunctionForArray(nicheFunction=conquadraticSkewed1, 
  #                      array=array(data= 1:10, dim =10), 
  #                      args=list(Xmin=0, Xmax=10, Xopt=5, Yopt=1))
}


nicheFunctionForRasterLayer <- function(nicheFunction, rasterLayer, args){ 
  # Function to apply a niche function over a raster layer
  #
  # Args:
  #   nicheFunction: the name of the niche function which is called
  #   rasterLayer: the rasterLayer of environmental values used to compute niche function.
  #   args : a list of the arguments of the niche function
  # 
  # Returns:
  #   A raster with values corresponding to the norm reaction
  values(rasterLayer) <- apply(X=as.array(getValues(rasterLayer)), MARGIN=1, FUN=nicheFunctionForValue, 
                      nicheFunction=nicheFunction, 
                      args=args)
  return(rasterLayer)
  # Ex :
  # raster <- raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)
  # nicheFunctionForRasterLayer(nicheFunction=conquadraticSkewed1, 
  #                            rasterLayer=raster, 
  #                            args=list(Xmin=0, Xmax=10, Xopt=5, Yopt=1))
  
}

geometricMean = function(x, na.rm=TRUE){
  #  A vectorized, zero- and NA-tolerant function for calculating geometric mean in R. The verbose mean calculation involving length(x) is necessary for the cases where x contains non-positive values.
  #
  # Args:
  #   x: an array
  #
  # Returns:
  #   the geometric mean of the array
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

nicheFunctionForRasterStack <- function(functionList, rasterStack, args){
  # Function to apply various niche functions to each layer of a rasterStack, combining them with a geometric mean. 
  #
  # Args:
  #   nichesFunctions: an ordered list containing the names of the niche functions which are called. Must be in the same order as rasterLayer and args
  #   rasterStack: the rasterstack of environmental values used to compute niche function.
  #   args : an ordered list containing  the lists of arguments necessary to call each niche function. 
  # 
  # Returns:
  #   A raster with values corresponding to the norm reaction
  
  # rearrange the list to get each element to be applied in an apply in the same index level
  X <- lapply(X=1:length(names(rasterStack)), 
              FUN=function(i, rasterStack){assign(paste("r",i), list(rasterStack[[i]], functionList[[i]], args[[i]])) },
              rasterStack)
  
  # Apply over each layer
  reactionNorm <- lapply(X, function(x){do.call(x[[2]], c(x[[1]], x[[3]]))})
  
  # Transform rasters to matrix
  unroll <- sapply(X=reactionNorm, FUN=getValues)
  
  # Combine the response with a geometric mean : gives a vector
  response <- apply(X=unroll, FUN=geometricMean, MARGIN=1)
  
  # Format the output
  rasterModel <- rasterStack[[1]]
  names(rasterModel) <- "Response"
  values(rasterModel) <- response

  return(rasterModel)
  
  # Ex:*
  # functionList <- list(conquadraticSkewed1, linearPositiveTwoParameters)
  # Data2 <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 
  # rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1),
  #                                         "BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)))
  # args <- list(list(Xmin=0, Xmax=10, Xopt=5, Yopt=1),list(X0=0, slope=1/2))
  # nicheFunctionForRasterStack(functionList, rasterStack, args)
  
}

nicheFunctionForArray2 <- function(functionList, Array, args){
  # Function to apply various niche functions to each column of an array, combining them with a geometric mean. 
  #
  # Args:
  #   nichesFunctions: an ordered list containing the names of the niche functions which are called. Must be in the same order as rasterLayer and args
  #   Array: the array of environmental values used to compute niche function.
  #   args : an ordered list containing  the lists of arguments necessary to call each niche function. 
  # 
  # Returns:
  #   A raster with values corresponding to the norm reaction
  
  # rearrange the list to get each element to be applied in an apply in the same index level
  X <- lapply(X=1:ncol(Array), 
              FUN=function(i, Array){assign(paste("r",i), list(Array[,i], functionList[[i]], args[[i]])) },
              Array)
  
  # Apply over each column
  reactionNorm <- lapply(X, function(x){do.call(x[[2]], c(x[[1]], x[[3]]))})
  
  # Transform rasters to matrix
  unroll <- sapply(X=reactionNorm, FUN=getValues)
  
  # Combine the response with a geometric mean : gives a vector
  response <- apply(X=unroll, FUN=geometricMean, MARGIN=1)
  return(response)
  
  # Ex:*
  # functionList <- list(conquadraticSkewed1, linearPositiveTwoParameters)
  # Data <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 
  # 
  # args <- list(list(Xmin=0, Xmax=10, Xopt=5, Yopt=1),list(X0=0, slope=1/2))
  # nicheFunctionForArray(functionList, Data, args)
  
}

################# End of Applying Niche Functions to Objects <<<<<<<<<<<<<<<<<<<<<<<

