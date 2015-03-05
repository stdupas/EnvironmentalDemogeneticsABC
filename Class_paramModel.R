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
    } else if(object@type_prior){
      return (TRUE)
    }
  }
  
)


#Initiateur with 1 arguments
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



#UserFriendly constructor with 1 arguments
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


findFunctionFromFile = function(model_type,fct_name){
  data_fct = read.table("functions.txt", sep = ";", header = TRUE)
  if(!is.element(fct_name, data_fct[,2])){
    stop("The function: ", fct_name," does not match any knowm function")
  } else {
    indice = match(fct_name, data_fct[,2])
    if(data_fct[indice,1]!=model_type){
      stop("This function is not a ", model_type, " function")
    } else {
      paul=unlist(strsplit((toString(data_fct[indice,4])),","))
      return (c(data_fct[indice,3], as.vector(paul)))
    }
  }
}




