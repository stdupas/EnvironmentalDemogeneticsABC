#Class of paramModel
setClass(
  Class = "ParamModel",
  representation = representation(
    name = "character",
    type_prior = "character",
    param_prior = "numeric"
    ),
  prototype = prototype(
    name = character(0),
    type_prior = character(0), 
    param_prior = numeric(0)
    ),
  validity = function(object){
    cat("---------- ParamModel : verification ----------\n")
    #Will be changed in order to check if type_prior match a known prior function
  
    if((is.null(object@type_prior))){
      stop("[ ParamModel : verificatiteston ] type_prior does not match any know prior function")
    } else {
      return (TRUE)
    }
  }
  
)

#Initiateur with 2 arguments
#setMethod(
#  f = "initialize",
#  signature = "ParamModel",
#  definition = function(.Object, type_prior, param_prior){
#    cat("---------- ParamModel : initiation ----------\n")
#    if(!missing(type_prior) && !(missing(param_prior))){
#      .Object@type_prior = type_prior
#      .Object@param_prior = param_prior
#      validObject(.Object)
#    } else {
#      stop("[ ParamModel : initiation ] Argument is missing")
#    }
#    return(.Object)
#  }
#)


#Initiateur with 0 arguments
setMethod(
  f = "initialize",
  signature = "ParamModel",
  definition = function(.Object, model_num){
    cat("---------- ParamModel : initiation ----------\n")
    .Object@name=c("model:",model_num)
    validObject(.Object)
    return(.Object)
  }
)



#UserFriendly constructor with 2 arguments
#paramModel = function(type_prior, param_prior){
#  cat("---------- ParamModel : construction ----------\n")
#  new(Class = "ParamModel", type_prior=type_prior, param_prior=param_prior)
#}

#UserFriendly constructor with 0 arguments
paramModel = function(model_num){
  cat("---------- ParamModel : construction ----------\n")
  new(Class = "ParamModel", model_num = model_num)
}

#Function to get the "type_prior" attribut
setGeneric("getType_prior",
           function(object){standardGeneric("getType_prior")})

setMethod("getType_prior", "ParamModel",
          function(object){
            return(object@type_prior)
          }
)

#Function to get the "param_prior" attribut
setGeneric("getParam_prior",
           function(object){standardGeneric("getParam_prior")})

setMethod("getParam_prior", "ParamModel",
          function(object){
            return(object@param_prior)
          }
)