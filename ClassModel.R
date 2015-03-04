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
    type_model = character(0), 
    param_model = NULL
  ),
)