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
)