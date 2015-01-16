
listOfNicheParameters <- function(rasterStack){
  # Ask to user for niche function to be applied to environemental variables, then append a list according to their parameters
  #
  # Args: 
  #   rasterStack: a rasterStack of environmental variables over wich we apply different niche functions for each layer.
  #
  # Returns:
  #   A list representing the decision tree for parameters
  
  Niche <- list()
  for(envrt_variable in names(rasterStack)) # envrt_variable <- "BIO1"
  {
    # add a box "BIO1" in simulation$Niche
    Niche[[envrt_variable]] <- list()
    
    # Ask to the user which niche model to apply (it has to be a function name !!!)
    niche_model <- readline(paste("Which niche model do you want for",envrt_variable,"? ")) # TODO : conquadraticsquewed
    
    for(niche_parameter in names(formals(niche_model))) # niche_parameter="X"
    {
      # add a box "X in simulations$Niche$BIO1
      Niche[[envrt_variable]][[niche_parameter]] <- list()
      
      # Ask user for prior distribution (it has to be a function name !!!)
      prior_distribution <- readline(paste("Which prior model do you want for",niche_parameter,"? ")) # TODO : uniform
      
      for(prior_parameter in names(formals(prior_distribution))) # prior_parameter = "n"
      {
        # add a box "n" in simulation$Niche$BIO1$X
        Niche[[envrt_variable]][[niche_parameter]][[prior_parameter]] <- numeric()
        
      }# end of loop over priors parameters
    }# end of loops over niche function parameters
  }# end of loop over environmental variable
  return(Niche)
}

listOfDispersionParameters <- function(){
  # Ask to user for dispersion function to be applied, then append a list according to its parameters
  # 
  # Args: takes no arguments
  #
  # Returns : 
  #   A list representing the decision tree for parameters
  Dispersion <- list()
  # Ask for dispersion model :
  dispersion_model <- readline(paste("Which dispersion model do you want ?"))
  
  # Get the arguments of the function and loop over them ((it has to be a function name !!!))
  for(dispersion_parameter in names(formals(dispersion_model))){
    Dispersion[[dispersion_parameter]] <- list()
    
    # Ask for prior distribution 
    prior_distribution <- readline(paste("Which prior model do you want for",dispersion_parameter,"? ")) # TODO : uniform
    # and loop over the parameters of prior distribution
    for(prior_parameter in names(formals(prior_distribution))){
      Dispersion[[dispersion_parameter]][[prior_parameter]] <- numeric()
    }# end of loop over prior parameters
  }# end of loop over dispersion parameters
  return(Dispersion)
}

listOfMutationParameters <- function(){
  # Ask to user for mutation function to be applied, then append a list according to its parameters
  # 
  # Args: takes no arguments
  #
  # Returns : 
  # A list representing the decision tree for parameters.
  Mutation <- list()
  
  # Ask for mutation model :
  mutation_model <- readline(paste("Which mutation model do you want ? "))
  
  # Get the arguments of the mutation function and loop over them ((it has to be a function name !!!))
  for(mutation_parameter in names(formals(mutation_model))){
    Mutation[[mutation_parameter]] <- list()
    
    # Ask for prior model 
    prior_distribution <- readline(paste("Which prior model do you want for",mutation_parameter,"? ")) # TODO : uniform
    # and loop over the parameters of the prior distribution
    for(prior_parameter in names(formals(prior_distribution))){
      Mutation[[mutation_parameter]][[prior_parameter]] <- numeric()
    }# end of loop over prior parameters
    
  }# end of loop over mutation model parameters
  return(Mutation)
}


askListOfParameters <- function(rasterStack){
  # Ask to user for niche, dispersion and mutation functions to be applied, appending a list with the models parameters
  # 
  # Args:
  #   rasterStack: a rasterStack of environmental variables over wich we apply different niche functions for each layer.
  #
  # Returns : 
  # A list representing the decision tree for parameters.
  Niche <- listOfNicheParameters(rasterStack=rasterStack)
  Dispersion <- listOfDispersionParameters()
  Mutation <- listOfMutationParameters()
  return(c(Niche,Dispersion,Mutation))
}
