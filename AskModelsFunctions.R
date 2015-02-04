
listOfNicheParameters <- function(rasterStack, nb_simulations){
  # Ask to user for niche function to be applied to environemental variables, then append a list according to their parameters
  #
  # Args: 
  #   rasterStack: a rasterStack of environmental variables over wich we apply different niche functions for each layer.
  #   nb_simulations: a numeric, the number of samples in the prior distribution of parameters
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
    
    # add a box "model"
    Niche[[envrt_variable]][["nicheModel"]] <- niche_model
    
    for(niche_parameter in names(formals(niche_model))[-1]) # niche_parameter="X"
    {
      # add a box "X in simulations$Niche$BIO1
      Niche[[envrt_variable]][[niche_parameter]] <- list()
      
      # Ask user for prior distribution (it has to be a function name !!!)
      prior_distribution <- readline(paste("Which prior model do you want for",niche_parameter,"? ")) # TODO : uniform
      
      # add a box "priorLaw"
      Niche[[envrt_variable]][[niche_parameter]][["priorLaw"]] <- prior_distribution
      
      for(prior_parameter in names(formals(prior_distribution))[-1]) # prior_parameter = "n"
      {
        # Ask user the values of prior parameters he wants
        parameter_value <- readline(paste("Which value for prior parameter", prior_parameter, "do you want ? ", sep=" "))
        
        # add a box "n" 
        Niche[[envrt_variable]][[niche_parameter]][[prior_parameter]] <- as.numeric(parameter_value)
        
      }# end of loop over priors parameters
      
      # add a box "values" and fill it taking the function and its arguments in the list
      Niche[[envrt_variable ]][[niche_parameter]][["Values"]] <- 
        do.call(what = Niche[[envrt_variable]][[niche_parameter]][["priorLaw"]],
                args = c(list(n = nb_simulations), Niche[[envrt_variable]][[niche_parameter]][-1]))
      
    }# end of loops over niche function parameters
  }# end of loop over environmental variable
  return(Niche)
}

listOfDispersionParameters <- function(nb_simulations){
  # Ask to user for dispersion function to be applied, then append a list according to its parameters
  # 
  # Args:
  #   nb_simulations: a numeric, the number of samples in the prior distribution of parameters
  #
  # Returns : 
  #   A list representing the decision tree for parameters
  
  Dispersion <- list()
  
  # Ask for dispersion model :
  dispersion_model <- readline(paste("Which dispersion model do you want ? "))
  # add a box "model"
  Dispersion[["DispersionModel"]] <- dispersion_model
  
  # Get the arguments of the function and loop over them ((it has to be a function name !!!))
  for(dispersion_parameter in names(formals(dispersion_model))[-1]){
    # add a box for the parameter
    Dispersion[[dispersion_parameter]] <- list()
    # Ask for prior distribution 
    prior_distribution <- readline(paste("Which prior model do you want for",dispersion_parameter,"? ")) # TODO : uniform
    # add a box "priorLaw"
    Dispersion[[dispersion_parameter]][["priorLaw"]] <- prior_distribution
    
    # and loop over the parameters of prior distribution
    for(prior_parameter in names(formals(prior_distribution))[-1]){
      # Ask user the values of prior parameters he wants and use it to fill the box
      parameter_value <- readline(paste("Which value for prior parameter", prior_parameter, "do you want ? ", sep=" "))
      Dispersion[[dispersion_parameter]][[prior_parameter]] <- as.numeric(parameter_value)
      
    }# end of loop over prior parameters
    
    # add a box "values" and fill it taking the function and its arguments in the list
    Dispersion[[dispersion_parameter]][["Values"]] <- 
      do.call(what = Dispersion[[dispersion_parameter]][["priorLaw"]],
              args = c(list(n = nb_simulations), Dispersion[[dispersion_parameter]][-1]))
    
  }# end of loop over dispersion parameters
  return(Dispersion)
}



listOfMutationParameters <- function(nb_simulations){ # nb_simulations <- 10
  # Ask to user for mutation function to be applied, then append a list according to its parameters
  # 
  # Args:
  #   nb_simulations: a numeric, the number of samples in the prior distribution of parameters
  #
  # Returns : 
  #    A list representing the decision tree for parameters.
  Mutation <- list()
  
  ##### Ask for mutation rate
  mutationRate <- list()
  mutationRatePrior <-  readline("Which prior for mutation rate do you want ? ")
  mutationRate[["priorLaw"]] <- mutationRatePrior
  for(prior_parameter in names(formals(mutationRate[["priorLaw"]]))[-1]){ # prior_parameter <- "min" prior_parameter <- "max"
    # Ask user the values of prior parameters he wants and use it to fill the box
    parameter_value <- readline(paste("Which value for prior parameter", prior_parameter, "do you want ? ", sep=" "))
    mutationRate[[prior_parameter]] <- as.numeric(parameter_value)
  }# end of loop over prior parameters
  mutationRate[["Values"]] <- do.call(what = mutationRate[["priorLaw"]],
            args = c(list(n = nb_simulations), mutationRate[-1]))
  
  ##### Ask for mutation model :
  mutation_model <- readline(paste("Which mutation model do you want ? "))
  # add a box for mutation model
  Mutation[["mutationModel"]] <- mutation_model
  
  # Get the arguments of the mutation function and loop over them ((it has to be a function name !!!))
  for(mutation_parameter in names(formals(mutation_model))[-1]){
    Mutation[[mutation_parameter]] <- list()
    
    # Ask for prior model 
    prior_distribution <- readline(paste("Which prior model do you want for",mutation_parameter,"? ")) # TODO : uniform
    # add a box for prior model
    Mutation[[mutation_parameter]][["priorLaw"]] <- prior_distribution
    
    # and loop over the parameters of the prior distribution
    for(prior_parameter in names(formals(prior_distribution))[-1]){
      # Ask user the values of prior parameters he wants and use it to fill the box
      parameter_value <- readline(paste("Which value for prior parameter", prior_parameter, "do you want ? ", sep=" "))
      Mutation[[mutation_parameter]][[prior_parameter]] <- as.numeric(parameter_value)
      
    }# end of loop over prior parameters
    
    Mutation[[mutation_parameter]][["Values"]] <- 
      do.call(what = Mutation[[mutation_parameter]][["priorLaw"]],
              args = c(list(n = nb_simulations), Mutation[[mutation_parameter]][-1]))
    
  }# end of loop over mutation model parameters
  
  # pass the sub-list to Mutation
  Mutation[["mutationRate"]] <- mutationRate
  return(Mutation)
}



askListOfParameters <- function(rasterStack, nb_simulations){
  # Ask to user for niche, dispersion and mutation functions to be applied, appending a list with the models parameters
  # 
  # Args:
  #   rasterStack: a rasterStack of environmental variables over wich we apply different niche functions for each layer.
  #
  # Returns : 
  # A list representing the decision tree for parameters.
  print("Please answer the questions for the carrying capacity niche model")
  NicheK <- listOfNicheParameters(rasterStack=rasterStack, nb_simulations = nb_simulations)
  print("Please answer the questions for the growth rate niche model")
  NicheR <- listOfNicheParameters(rasterStack=rasterStack, nb_simulations = nb_simulations)
  Dispersion <- listOfDispersionParameters(nb_simulations = nb_simulations)
  Mutation <- listOfMutationParameters(nb_simulations = nb_simulations)
  return(list("NicheK" = NicheK, "NicheR" = NicheR, "Dispersion" = Dispersion, "Mutation" = Mutation))
}

referenceTableFromList <- function(ParamList){
  # Makes a reference table for abc inference, with all values of parameters sampled from their prior distribution.
  # 
  # Args :
  #   ParamList : the list constructed when asking user with askListOfParameters function
  #
  # Returns :
  #   The reference table
  
  ### Find niche informations and prior sampling values :
  niche <- c()
  nicheNames <- c()
  
  # Loop over the layers of rasterStack
  for(bioLayer in names(ParamList[["Niche"]])){
    foo <- ParamList[["Niche"]][[bioLayer]][1]
    # Loop over the parameters of the model, omitting the first case (containing the name of the model used)
    for(parameter in names(ParamList[["Niche"]][[bioLayer]][-1])){
      # construct the full name of the parameter
      toto <- paste("Niche",bioLayer, foo, parameter, sep=".")
      
      niche <- rbind(niche, ParamList[["Niche"]][[bioLayer]][[parameter]][["Values"]])
      nicheNames <- c(nicheNames, toto)
    }
  }
  rownames(niche) <- nicheNames
  
  ### Find dispersion informations and prior sampling values :
  dispersion <- c()
  dispersionNames <- c()
  
  # Loop over the parameters of the model, omitting the first case (containing the name of the model used)
  for(parameter in names(ParamList[["Dispersion"]][-1])){
    # construct the full name of the parameter
    toto <- paste("Dispersion", ParamList[["Dispersion"]][[1]], parameter, sep=".")
    
    dispersion <- rbind(dispersion, ParamList[["Dispersion"]][[parameter]][["Values"]])
    dispersionNames <- c(dispersionNames, toto)
  }
  rownames(dispersion) <- dispersionNames
  
  ### Find mutation informations and prior sampling values :
  mutation <- c()
  mutationNames <- c()
  
  # Loop over the parameters of the model, omitting the first case (containing the name of the model used)
  for(parameter in names(ParamList[["Mutation"]][-1])){
    # construct the full name of the parameter
    toto <- paste("Mutation", ParamList[["Mutation"]][[1]], parameter, sep=".")
    
    mutation <- rbind(mutation, ParamList[["Mutation"]][[parameter]][["Values"]])
    mutationNames <- c(mutationNames, toto)
  }
  rownames(mutation) <- mutationNames
  
  ### Concatenate all
  return(rbind(niche, dispersion, mutation))
}


getFunctionListNiche <- function(ParamList, sublist){
  # Get a list of niche models functions stored in ParamList
  #
  # Args:
  #   ParamList: the list of model parameters created using askListOfParameters function
  #   sub: caracetr vector, the sublist (NicheR or NicheK) to explore
  #
  # Returns:
  #   the list of niche model functions used
  lapply(X = names(ParamList[[sublist]]), FUN = function(x, ParamList){ParamList[[sublist]][[x]][[1]]}, ParamList=ParamList)
}

getFunctionDispersion <- function(ParamList){
  # Get a list of dispersion models functions stored in ParamList
  #
  # Args:
  #   ParamList: the list of model parameters created using askListOfParameters function
  #
  # Returns:
  #   the list of niche model functions used
  return(ParamList[["Dispersion"]][[1]])
}


getArgsListNiche <- function(simulation, ParamList, sublist){
  # Get a list of parameters of niche models functions and their values stored in ParamList
  #
  # Args:
  #   simulations : the index in which to get the value stored in the vector Values (by prior sampling) in ParamList.
  #   ParamList: the list of model parameters created using askListOfParameters function
  #
  # Returns:
  #   the list of niche model functions parameters 
  
  # initialize the list of parameters
  argsList <- list()
  
  # Loop over the layers of environmental variables
  for(layer in names(ParamList[[sublist]])){
    # initialize the vectors
    argsValues <- c()
    argsNames <- c()
    
    # Loop over the parameters (omitting the first element of the list, containing the name of the model)
    for(param in names(ParamList[[sublist]][[layer]])[-1]){
      argsNames <- c(argsNames, as.character(param))
      argsValues <- c(argsValues, value=ParamList[[sublist]][[layer]][[param]][["Values"]][simulation])
    }
    names(argsValues)<- argsNames
    argsList[[length(argsList)+1]] <- as.list(argsValues)
  }
  return(argsList)
}

getArgsListDispersion <- function(simulation, ParamList){
  # Get a list of parameters of dispersion models functions and their values stored in ParamList
  #
  # Args:
  #   simulations : the index in which to get the value stored in the vector Values (by prior sampling) in ParamList.
  #   ParamList: the list of model parameters created using askListOfParameters function
  #
  # Returns:
  #   the list of niche model functions parameters 
  
  # initialize the list of parameters
  argsList <- list()
  
  # initialize the vectors
  argsValues <- c()
  argsNames <- c()
  
  # Loop over the parameters (omitting the first element of the list, containing the name of the model)
  for(param in names(ParamList[["Dispersion"]])[-1]){
    argsNames <- c(argsNames, as.character(param))
    argsValues <- c(argsValues, value=ParamList[["Dispersion"]][[param]][["Values"]][simulation])
  }
  names(argsValues)<- argsNames
  return(as.list(argsValues))
}