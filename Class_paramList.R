# Sources for needed classes
source("Class_composante.R")

# Class paramList
setClass(
    Class="ParamList", 
    representation=representation(
        niche_r="Composante",
        niche_k="Composante",
        dispersion="Composante",
        mutation="Composante",
        names_stacks = "character"
    )
)

#Function to save the ParamList
setGeneric("saveParamList",
           function(object, file){standardGeneric("saveParamList")})

setMethod("saveParamList", "ParamList",
          function(object, file){
              # Function to save a paramList object into a RDS file
              # Args:
              #     file: a string containing the name of the file where the paramList will be saved
              # Example:
              #     saveParamList("debug.rds")
              saveRDS(object, file)
          }
)


#Function to load the ParamList
loadParamList = function(file){
    # Function to load a paramList object from a RDS file
    # Args:
    #     file: a string containing the name of the file where the paramList will be load
    #Return:
    #   object: the paramList in the file
    # Example:
    #     loadParamList("debug.rds")
    object = readRDS(file)
    return(object)
}

#Function to get the "result_prior"
getResultPrior = function(object,comp,param){
    # Function that asks the user which prior results he wants and then returns it
    # Args:
    #     object: a paramList object where the wanted prior results are
    #     comp: the name of the composante where the wanted prior results are
    #     param: a vector containing the model number and the number of the parameter model
    #            where the wanted prior results are. or 0 if the wanted prior results are in
    #            the independant model
    #Return:
    #   object: a list contaning te prior results asked
    # Example:
    #     getResultPrior(my_paramlist, "niche_k", c(1,2))
    #       returns the prior results of the sencond hyper_parameter of the first model of the component niche_k 
    #       in the paramList my_paramlist
    
    if(object@method == "Bayesian" || object@method == "Likelihood") {
        cat("ERROR: Your method doesn't allow any prior result")
    } else if(comp == "niche_r") {
        if(param[1] == 0) {
            res = getResult_prior(object@niche_r@independance)
        } else {
            res = getResult_prior(object@niche_r@listModel[[param[1]]]@param_model[[param[2]]])
        }
    } else if(comp == "niche_k") {
        if(param[1] == 0) {
            res = getResult_prior(object@niche_k@independance)
        } else {
            res = getResult_prior(object@niche_k@listModel[[param[1]]]@param_model[[param[2]]])
        }
    } else if(comp == "dispersion") {
        res = getResult_prior(object@dispersion@listModel[[param[1]]]@param_model[[param[2]]])
    } else if(comp == "mutation") {
        res = getResult_prior(object@mutation@listModel[[param[1]]]@param_model[[param[2]]])
    } else if(comp == "generation") {
        if(param[1] == 0) {
            res = getResult_prior(object@generation@independance)
        } else {
            res = getResult_prior(object@generation@listModel[[param[1]]]@param_model[[param[2]]])
        }
    }
    if(length(res) == 0) {
        cat("ERROR: You did not compute the result prior yet.")
    } else {
        return(res)
    }
}


# Constructor of paramList
paramList = function(env) {
    # Function to initiate the construction of an object of classe paramList
    # Args:
    #     env: A rasterStack or an Array that contains the differnets environnement variables
    #Return:
    #   .Object: an object of class paramList   

    flag = -1
    while(flag == -1) {
        cat("[Type 0 to quit] What do you want to do?\n1: Forward\n2: Backward\n")
        scanner = as.integer(readline())
        if(is.na(scanner) || scanner>2 || scanner<0) {
            print("ERROR: Your entry is incorrect, please try again")
        } else if(scanner == 0) {
            stop("You have stopped the program")
        } else {
            flag = 1
        }
    }
    # if the environmental data is an array
    if(is.array(env)) {
        na = dimnames(env)[[3]]
    }
    # if the environmental data is a rasterstack
    else {
        na = names(env)
    }
    # initialisation of the object
    if(scanner == 1) {
        .Object = forward(na)
    } else if(scanner == 2) {
        .Object = backward(na)
    }
    return(.Object)
}

# Function to assess the priors
setResultPrior = function(object) {
    # Function that asks the user which prior results he wants to generate
    # Args:
    #     object: a paramList object
    #Return:
    #   setResultPrior...(object): a paramList containing the same informations that object
    #                              but also a full result_prior slot
    #EXample:
    #   my_paramlist = setResultPrior(my_paramlist)
    `
    if(class(object)[1] == "Backward") {
        return(setResultPriorBack(object))
    } else if(class(object)[1] == "Forward") {
        return(setResultPriorFor(object))
    } else {
        stop("The specified object is not valid.")
    }
}

# Function to change the parameters in paramList
setParamList = function(object) {
    # Function that asks the user whhat he wants to change
    # Args:
    #     object: a paramList object
    #Return:
    #   set...(object): a paramList which is the modified object
    #EXample:
    #   my_paramlist = setParamList(my_paramlist)
    
    if(class(object)[1] == "Backward") {
        return(setBackward(object))
    } else if(class(object)[1] == "Forward") {
        return(setForward(object))
    } else {
        stop("The specified object is not valid.")
    }
}
