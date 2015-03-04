#Class of Model

source("ClassParamModel.R")

setClass(
  #exemple of creation: 
  #new(Class="Model", type_model="paul", param_model = new(Class="ParamModel", type_prior="julie", param_prior=c(0,1,2,3)))
  Class = "Model",
  representation = representation(
    type_model = "character",
    param_model = "list"
  ),
  prototype = prototype(
    type_model = "vide", 
    param_model = NULL
  )
  #validity = function(object){ ## object : type_model doesn't match any known function
  #if (object@type_model == "vide" || autre types de fonctions){
  #  return(FALSE)
  #} else {
  #  return(TRUE)
  #}
  #}
)