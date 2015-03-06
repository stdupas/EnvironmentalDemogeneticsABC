#Class of Model

source("Class_paramModel.R")

setClass(
  Class = "Model",
  representation = representation(
    name = "character",
    type_model = "character",
    param_model = "list"
  ),
  prototype = prototype(
    name = character(0),
    type_model = character(0), 
    param_model = list(0)
  ),
  validity = function(object){
    #Will be changed in order to check if type_model match a known model function
    if(is.null(object@type_model)){
      stop("[ Model : verification ] type_model does not match any know model function")
    } else{
    
      return (TRUE)
    }
  }
  
)

#Initiateur
setMethod(
  f = "initialize",
  signature = "Model",
  definition = function(.Object, composante_name, model_num){
    .Object@name=c(composante_name, model_num)
    print("[Type 0 to exit] What function do you want to use for the model ?")
    data_fct = read.table("functions.txt", sep = ";", header = TRUE, as.is=rep(TRUE, 4))
    possible = which(data_fct[,1]==composante_name)

    # print the possible functions in order (with a number to know which one to use)
    possible_number = 1:length(possible)
    possible_function = data_fct[possible,2]
    possible_print = paste(possible_number, rep(":",length(possible)), possible_function)
    print(possible_print)

    # choose the function with the given number
    # repeat while the given number is incorrect
    flag = -1
    while (flag == -1) {
        choice_number = as.integer(readline())
        if(choice_number!=0 && choice_number<=length(possible) && choice_number>0 && !is.na(choice_number)) {
            choice_function = possible_function[choice_number]
            .Object@type_model = toString(choice_function)
            flag = 1
        }
        else if(choice_number==0 && !is.na(choice_number)) {
            stop("Stop the program.")
        }
        else {
            print("Wrong number, please type a number in the list :")
            print(possible_print)
        } 
    }

    # ask which parameters and prior functions they want
    vec = findFunctionFromFile(composante_name, .Object@type_model)
    mod = NULL
    for(i in 1:vec[1]){
      print(paste("========== ParamModel : ",.Object@type_model,", parameter: ",vec[i+1]," =========="))
      mod = c(mod, paramModel(model_num))
    }
    .Object@param_model = mod  
    validObject(.Object)
    return(.Object)
  }
)


#UserFriendly constructor
#model = function(type_model, param_model){
#  cat("---------- Model : construction ----------\n")
#  new(Class = "Model", type_model=type_model, param_model=param_model)
#}

model = function(composante_name, model_num){
  new(Class = "Model", composante_name=composante_name, model_num=model_num)
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

# Function to update the number of models
setGeneric(
    name="setNumModel",
    def=function(object, nb) {standardGeneric("setNumModel")}
)


setMethod(
    f="setNumModel",
    signature="Model",
    definition=function(object, nb) {
        object@name[2]=nb
        return(object)
    }
)


