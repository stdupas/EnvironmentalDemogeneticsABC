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
              saveRDS(object, file)
          }
)


#Function to load the ParamList
loadParamList = function(file){
    object = readRDS(file)
    return(object)
}

#Function to get the "result_prior"
#Args: object paramList, component name, c(model number, parameter model) or 0 for independant model
setGeneric("getResultPrior",
           function(object,comp,param){standardGeneric("getResultPrior")})

setMethod("getResultPrior", "ParamList",
          function(object,comp,param){
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
)

# Constructor of paramList
paramList = function(env) {
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

