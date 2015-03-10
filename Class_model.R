#Class of Model

source("Class_paramModel.R")

setClass(
  Class = "Model",
  representation = representation(
      method = "character",
    name = "character",
    type_model = "character",
    param_model = "list",
    param_name = "character"
  ),
  prototype = prototype(
    name = character(0),
    type_model = character(0), 
    param_model = list(0),
    param_name = character(0)
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
  definition = function(.Object, composante_name, model_num, method){
    .Object@name=c(composante_name, model_num)
    .Object@method = method
    print("[Type 0 to exit] What function do you want to use for the model ?")
    data_fct = read.table("functions.txt", sep = ";", header = TRUE, as.is=rep(TRUE, 4))
    if(composante_name == "niche_r" || composante_name == "niche_k") {
        possible = which(data_fct[,1]=="niche")
    } else {
        possible = which(data_fct[,1]==composante_name)
    }
    

    # print the possible functions in order (with a number to know which one to use)
    possible_number = 1:length(possible)
    possible_function = data_fct[possible,2]
    possible_print = paste(possible_number, rep(":",length(possible)), possible_function)

    # choose the function with the given number
    # repeat while the given number is incorrect
    flag = -1
    while (flag == -1) {
      print(possible_print)
      choice_number = as.integer(readline())
      if(choice_number!=0 && choice_number<=length(possible) && choice_number>0 && !is.na(choice_number)) {
        choice_function = possible_function[choice_number]
        .Object@type_model = toString(choice_function)
        flag = 1
      }
      else if(choice_number==0 && !is.na(choice_number)) {
        stop("You have stopped the program")        }
      else {
        print("ERROR: Your entry is incorrect, please try again")
      } 
    }

    # ask which parameters and prior functions they want
    if(composante_name == "niche_r" || composante_name == "niche_k") {
        vec = findFunctionFromFile("niche", .Object@type_model)
    } else {
        vec = findFunctionFromFile(composante_name, .Object@type_model)
    }
    
    mod = NULL
    param_name = NULL
    for(i in 1:vec[1]){
      print(paste("========== ParamModel : ",.Object@type_model,", parameter: ",vec[i+1]," =========="))
      mod = c(mod, paramModel(model_num, .Object@method))
      param_name = c(param_name, vec[i+1])
    }
    .Object@param_model = mod 
    .Object@param_name = param_name 
    validObject(.Object)
    return(.Object)
  }
)


model = function(composante_name, model_num, method){
  new(Class = "Model", composante_name=composante_name, model_num=model_num, method=method)
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

# Function to print the parameters of the model functions
setMethod(
    f="show", 
    signature="Model",
    definition=function(object) {
        for(i in 1:length(object@param_name)) {
            cat("          ",object@param_name[i],": ")
            print(object@param_model[[i]])
        }
    }
)

# Function to change the type of model (function to use for the model)
setGeneric(
    name="setTypeModel",
    def=function(object) {standardGeneric("setTypeModel")}
)


setMethod(
    f="setTypeModel",
    signature="Model",
    definition=function(object) {
        cat("The actual type of function for this model is :", object@type_model,"\n")
        # Ask what type the user wants to use
        newObject = model(object@name[1],object@name[2],object@method)
        return(newObject)
    }
)

# Function to change a prior
setGeneric(
    name="setPrior",
    def=function(object) {standardGeneric("setPrior")}
)


setMethod(
    f="setPrior",
    signature="Model",
    definition=function(object) {
        flag = -1
        while(flag == -1) {
          # Ask which parameter the user wants to change
          cat("[Type 0 to quit] Which parameter do you want to change ?\n")
          cat(paste(1:length(object@param_name),":",object@param_name))
            choice_number = as.integer(readline())
            if(choice_number!=0 && choice_number<=length(object@param_name) && choice_number>0 && !is.na(choice_number)) {
              flag2 = -1
              while(flag2 == -1) {
                # Ask if the user wants to change the prior function or only a parameter of a prior function
                cat("[Type 0 to quit] Change the prior function or the value of a parameter value ?\n")
                cat("1: Prior function\n2: Parameter value")
                choice_number2 = as.integer(readline())
                if(choice_number2!=0 && choice_number2<=2 && choice_number2>0 && !is.na(choice_number2)) {
                  if(choice_number2==1) {
                    object@param_model[[choice_number]] = setType_prior(object@param_model[[choice_number]])
                  }
                  else if(choice_number2==2) {
                    object@param_model[[choice_number]] = setParam_prior(object@param_model[[choice_number]])
                  }                    
                  flag2 = 1
                }
                else if(choice_number2==0 && !is.na(choice_number2)) {
                  stop("You have stopped the program")
                }else {
                  print("ERROR: Your entry is incorrect, please try again")
                }
              }
              flag = 1
            }
            else if(choice_number==0 && !is.na(choice_number)) {
              stop("You have stopped the program")
            }
            else {
              print("ERROR: Your entry is incorrect, please try again")
            }
        }
        return(object)
    }
)

# Function to assess the prior values for one model
setGeneric(
    name="setResultPriorMod",
    def=function(object, all) {standardGeneric("setResultPriorMod")}
)


setMethod(
    f="setResultPriorMod",
    signature="Model",
    definition=function(object, all) {
        if(all == 0) {
            flag = -1
            while(flag == -1) {
              # ask which parameter the user wants to assess the prior values
              cat("[Type 0 to quit] Which parameter do you want to assess? \n")
              cat(paste(1:length(object@param_name),":",object@param_name,"\n"))
              cat(length(object@param_name)+1,":", "all parameters\n")
              choice = as.integer(readline())
              if(choice!=0 && choice<=length(object@param_name)+1 && choice>0 && !is.na(choice)) {
                if(choice == length(object@param_name)+1) {
                  for(i in 1:length(object@param_name)) {
                    object@param_model[[i]] = setResult_prior(object@param_model[[i]])
                  }
                }
                else {
                  object@param_model[[choice]] = setResult_prior(object@param_model[[choice]])                  
                }
                flag = 1
              }
              else if (choice == 0 && !is.na(choice)) {
                stop("You have stopped the program")                }
              else {
                print("ERROR: Your entry is incorrect, please try again")
              }
            }
        } else if(all == 1) {
          for(i in 1:length(object@param_name)) {
            object@param_model[[i]] = setResult_prior(object@param_model[[i]])
          }
        }
        return(object)
    }
)