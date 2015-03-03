#Class of paramModel
setClass(
  Class = "ParamModel",
  representation = representation(
    type_prior = "character",
    param_prior = "numeric"
    ),
  prototype = prototype(
    type_prior = "vide", 
    param_prior = NULL
    )
  #validity = function(object){ ## object : type_prior doesn't match any known function
    #if (object@type_prior == "vide" || autre types de fonctions){
    #  return(FALSE)
    #} else {
    #  return(TRUE)
    #}
  #}
)