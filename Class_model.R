############### CLASS: "Model" ###############

source("Class_paramModel.R")

############### Creation of class "Model" ###############
    # method: a string containing the method used by the paramList (ABC, Bayesian, Likelihood)
    # name: a vector containing (component name, number of the model, model name based on names from stack or array)
    # type_model: a string containing the function used for the model
    # param_model: a list containing objects of class ParamModel (an object for each parameters of the model function)
    # param_name: a vector containing the names of the function's parameters
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
        # Check if type_model match a known model function
        if(is.null(getType_model(object))){
            stop("[ Model : verification ] type_model does not match any know model function")
        } else{
            
            return (TRUE)
        }
    }
    
)


############### Initiator ###############
    # Function to initialize an object of class "Model"
    # Args:
    #       composante_name: a string containing the name of the component where the model is 
    #           (niche_r, niche_k, dispersion, mutation, generation)
    #       model_num: an integer containing the number of the actual model
    #       method: a string containing the method used by the paramList (ABC, Bayesian, Likelihood)
    #       nameStack: name of the model from array or rasterStack, or "STANDARD"
    #           if the component is independant of environmental data
    # Example:
    #       myModel = model("niche_r", 1, "ABC", "BIO12")
setMethod(
    f = "initialize",
    signature = "Model",
    definition = function(.Object, composante_name, model_num, method, nameStack){
        .Object@name=c(composante_name, model_num, nameStack)
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
            vec = findFunctionFromFile("niche", getType_model(.Object))
        } else {
            vec = findFunctionFromFile(composante_name, getType_model(.Object))
        }
        
        mod = NULL
        param_name = NULL
        for(i in 1:vec[1]){
            print(paste("========== ParamModel : ",getType_model(.Object),", parameter: ",vec[i+1]," =========="))
            mod = c(mod, paramModel(model_num, getMethodeMod(.Object)))
            param_name = c(param_name, vec[i+1])
        }
        .Object@param_model = mod 
        .Object@param_name = param_name 
        validObject(.Object)
        return(.Object)
    }
)
    # User-friendly constructor
model = function(composante_name, model_num, method, nameStack){
    new(Class = "Model", composante_name=composante_name, model_num=model_num, method=method, nameStack = nameStack)
}


#############################################
############### Get functions ###############
#############################################

############### getType_Model ###############
    # Function to get the "type_model" attribut (model function's name)
setGeneric("getType_model",
           function(object){standardGeneric("getType_model")})

setMethod("getType_model", "Model",
          function(object){
              return(object@type_model)
          }
)

############### getParam_nameMod ###############
    # Function to get the "param_name" attribut (parameters' names)
setGeneric("getParam_nameMod",
           function(object){standardGeneric("getParam_nameMod")})

setMethod("getParam_nameMod", "Model",
          function(object){
              return(object@param_name)
          }
)

############### getNameModel ###############
    # Function to get the "name" attribut
setGeneric("getNameModel",
           function(object){standardGeneric("getNameModel")})

setMethod("getNameModel", "Model",
          function(object){
              return(object@name)
          }
)

############### getMethodMod ###############
    # Function to get the "method" attribut
setGeneric("getMethodeMod",
           function(object){standardGeneric("getMethodeMod")})

setMethod("getMethodeMod", "Model",
          function(object){
              return(object@method)
          }
)


#############################################
############### Set functions ###############
#############################################

############### setNumModel ###############
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

############### setTypeModel ###############
    # Function to change the type of model
setGeneric(
    name="setTypeModel",
    def=function(object) {standardGeneric("setTypeModel")}
)

setMethod(
    f="setTypeModel",
    signature="Model",
    definition=function(object) {
        cat("The actual type of function for this model is :", getType_model(object),"\n")
        # Ask what type the user wants to use
        newObject = model(getNameModel(object)[1],getNameModel(object)[2],getMethodeMod(object))
        return(newObject)
    }
)

############### setPrior ###############
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
            cat(paste(1:length(getParam_nameMod(object)),":",getParam_nameMod(object)))
            choice_number = as.integer(readline())
            if(choice_number!=0 && choice_number<=length(getParam_nameMod(object)) && choice_number>0 && !is.na(choice_number)) {
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

############### setResultPrior ###############
    # Function to compute the prior values for one model
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
                cat(paste(1:length(getParam_nameMod(object)),":",getParam_nameMod(object),"\n"))
                cat(length(getParam_nameMod(object))+1,":", "all parameters\n")
                choice = as.integer(readline())
                if(choice!=0 && choice<=length(getParam_nameMod(object))+1 && choice>0 && !is.na(choice)) {
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
            for(i in 1:length(getParam_nameMod(object))) {
                object@param_model[[i]] = setResult_prior(object@param_model[[i]])
            }
        }
        return(object)
    }
)


#############################################
############### Show functions ##############
#############################################

############### show ###############
    # Function to print the parameters of the model functions
setMethod(
    f="show", 
    signature="Model",
    definition=function(object) {
        for(i in 1:length(getParam_nameMod(object))) {
            cat("          ",getParam_nameMod(object)[i],": ")
            print(object@param_model[[i]])
        }
    }
)

