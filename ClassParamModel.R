#Class of paramModel
setClass(
  Class = "ParamModel",
  representation = representation(
    type_prior = "character",
    param_prior = "numeric"
    ),
  prototype = prototype(
    type_prior = character(0), 
    param_prior = numeric(0)
    ),
  validity = function(object){
    cat("---------- ParamModel : verification ----------\n")
    #Will be changed in order to check if type_prior match a known prior function
    if(object@type_prior == "Not a function"){
      stop("[ ParamModel : verification ] type_prior does not match any know prior function")
    } esle {
      return (TRUE)
    }
  }
  
)



#UserFriendly constructor
paramModel = function(type_prior, param_prior)


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