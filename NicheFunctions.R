################### Formating Niche Functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
  return(y)
}

################# End of formatted Niche Functions <<<<<<<<<<<<<<<<<<<<


################# Applying Niche Functions to Objects >>>>>>>>>>>>>>>>>>

nicheFunctionForValue <- function(nicheFunction, x, args){
  # Function to apply a niche function to a single value
  #
  # Args:
  #   nicheFunction: the name of the niche function which is called
  #   rasterStack: the rasterStack of environmental values used to compute niche function.
  #   args : a list of the arguments of the niche function
  # 
  # Returns:
  #   The value corresponding to the norm reaction
  args <- c(list(x), args)
  do.call(nicheFunction, args)
  #Ex : nicheFunctionForValue(conquadraticSkewed1, x=4, args=list(Xmin=0, Xmax=10, Xopt=5, Yopt=1))
}

nicheFunctionForArray <- function(nicheFunction, array, args){ 
  # Function to apply a niche function to an array
  #
  # Args:
  #   nicheFunction: the name of the niche function which is called
  #   array: the array of environmental values used to compute niche function.
  #   args : a list of the arguments of the niche function
  # 
  # Returns:
  #   An array corresponding to the norm reaction
  
  # raster <- raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)
  # nicheFunction <- conquadraticSkewed1
  
  apply(X=array, MARGIN=1, FUN=nicheFunctionForValue, nicheFunction=nicheFunction, args=args)
  # nicheFunctionForArray(nicheFunction=conquadraticSkewed1, 
  #                      array=array(data= 1:10, dim =10), 
  #                      args=list(Xmin=0, Xmax=10, Xopt=5, Yopt=1))
}


# Applying to rasterStacks
nicheFunctionForRaster <- function(nicheFunction, raster, ...){ 
  # Function to apply a niche function to a rasterStack.
  #
  # Args:
  #   nicheFunction: the name of the niche function which is called
  #   rasterStack: the rasterStack of environmental values used to compute niche function.
  #   ... : the arguments of the niche function
  # 
  # Returns:
  #   The raster corresponding to the norm reaction
  
  # raster <- raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)
  # nicheFunction <- conquadraticSkewed1
  
  apply(X=as.array(values(raster), MARGIN=1, FUN=nicheFunctionForValue, nicheFunction=nicheFunction, ...))
}

################# End of Applying Niche Functions to Objects <<<<<<<<<<<<<<<<<<<<<<<

# Non function stepwise model, has to be a bit modified cause od "coaltable" missing in arguments
stepwise <- function(initial_genetic_value,stepvalue)
{
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

#1: sigmaDisp    2: gammaDisp
fat_tail2 <- function(x,sigma,gamma){ x^gamma*exp(-2*x/(sigma^0.5)) }

# prior function to call to get its parameters (with uniform it's not really interesting, but with others...)
uniform <- function(n, min, max){runif(n=n, min=min,max=max)}


###### Create and append list
simulations <-list(Niche=list(), Dispersion=list(), Mutation=list())

# Append "Niche" part of the list
for(envrt_variable in names(rasterStack)) # envrt_variable <- "BIO1"
{
  # add a box "BIO1" in simulation$Niche
  simulations[["Niche"]][[envrt_variable]] <- list()
  
  # Ask to the user which niche model to apply (it has to be a function name !!!)
  niche_model <- readline(paste("Which niche model do you want for",envrt_variable,"? ")) # TODO : conquadraticsquewed
  
  for(niche_parameter in names(formals(niche_model))) # niche_parameter="X"
  {
    # add a box "X in simulations$Niche$BIO1
    simulations[["Niche"]][[envrt_variable]][[niche_parameter]] <- list()
    
    # Ask user for prior distribution (it has to be a function name !!!)
    prior_distribution <- readline(paste("Which prior model do you want for",niche_parameter,"? ")) # TODO : uniform
    
    for(prior_parameter in names(formals(prior_distribution))) # prior_parameter = "n"
    {
      # add a box "n" in simulation$Niche$BIO1$X
      simulations[["Niche"]][[envrt_variable]][[niche_parameter]][[prior_parameter]] <- numeric()
      
    }
  }
}

# Append "Dispersion" part of the list
dispersion_model <- readline(paste("Which dispersion model do you want ?"))
for(dispersion_parameter in names(formals(dispersion_model))){
  simulations[["Dispersion"]][[dispersion_parameter]] <- list()
  prior_distribution <- readline(paste("Which prior model do you want for",dispersion_parameter,"? ")) # TODO : uniform
  
  for(prior_parameter in names(formals(prior_distribution))){
    simulations[["Dispersion"]][[prior_parameter]] <- numeric()
  }
}

# Append "Mutation" part of the list
mutation_model <- readline(paste("Which mutation model do you want ? "))
for(mutation_parameter in names(formals(mutation_model))){
  simulations[["Mutation"]][[mutation_parameter]] <- list()
  prior_distribution <- readline(paste("Which prior model do you want for",mutation_parameter,"? ")) # TODO : uniform
  
  for(prior_parameter in names(formals(prior_distribution))){
    simulations[["Mutation"]][[prior_parameter]] <- numeric()
  }
}