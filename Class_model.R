#Class of Model

source("Class_paramModel.R")

setClass(
  Class = "Model",
  representation = representation(
    type_model = "character",
    param_model = "list"
  ),
  prototype = prototype(
    type_model = character(0), 
    param_model = NULL
  ),
  validity = function(object){
    cat("---------- Model : verification ----------\n")
    #Will be changed in order to check if type_model match a known model function
    if(object@type_model == "Not a function"){
      stop("[ Model : verification ] type_model does not match any know model function")
    } else {
      return (TRUE)
    }
  }
  
)

#Initiateur
setMethod(
  f = "initialize",
  signature = "Model",
  definition = function(.Object, type_model, param_model){
    cat("---------- Model : initiation ----------\n")
    if(!missing(type_model) && !(missing(param_model))){
      .Object@type_model = type_model
      .Object@param_model = param_model
      validObject(.Object)
    } else {
      stop("[ Model : initiation ] Argument is missing")
    }
    return(.Object)
  }
)


#UserFriendly constructor
model = function(type_model, param_model){
  cat("---------- Model : construction ----------\n")
  new(Class = "Model", type_model=type_model, param_model=param_model)
}


#Function to get the "type_model" attribut
setGeneric("getType_model",
           function(object){standardGeneric("getType_model")})

setMethod("getType_model", "Model",
          function(object){
            return(object@type_model)
          }
)

#Function to get the "param_model" attribut
setGeneric("getParam_model",
           function(object){standardGeneric("getParam_model")})

setMethod("getParam_model", "Model",
          function(object){
            return(object@param_model)
          }
)