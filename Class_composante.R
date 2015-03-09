# Sources for needed classes
source("Class_model.R")

# Class composante
setClass(
    Class="Composante", 
    representation=representation(
        name="character",
        listModel="list",
        nbModel="numeric"
    ),
    prototype=prototype(
        name=character(0),
        listModel=list(0),
        nbModel=numeric(0)
    ),
    validity=function(object) {
        if(object@nbModel <= 0) {
            stop("[Composante validation] Number of models is not positive\n")
        }
        if(length(object@listModel)!=object@nbModel) {
            stop("[Composante validation] Number of models is different from nbModel\n")
        }
        return(TRUE)
    }
)

# Constructor of composante
# setMethod(
#     f="initialize",
#     signature="Composante",
#     definition=function(.Object) {
#         cat("---------- Composantes : initiation ----------\n")
#         if(!missing(listModel) && !missing(nbModel)) {
#             .Object@listModel = listModel
#             .Object@nbModel = nbModel
#             validObject(.Object)
#         }
#         else {
#             stop("[Composante initiation] Missing argument(s)\n")
#         }
#         return(.Object)
#     }
# )

setMethod(
    f="initialize",
    signature="Composante",
    definition=function(.Object, name) {
        .Object@name=name
        
        # repeat while the given number is incorrect
        flag = -1
        while(flag == -1) {
            choice_number = as.numeric(readline(paste("How many models for",name,"? (press 0 to quit)")))
            if(choice_number!=0 && choice_number>0 && !is.na(choice_number)) {
                .Object@nbModel = choice_number
                flag = 1
            }
            else if(choice_number==0 && !is.na(choice_number)) {
                stop("Stop the program.")
            }
            else {
                print("Wrong number.")
            } 
        }

        mod = NULL
        for(i in 1:.Object@nbModel) {
            print(paste("========== Composante : ",name,", model n°", i," =========="))
            mod = c(mod, model(name,i))
        }

        .Object@listModel = mod
        validObject(.Object)
        return(.Object)
    }
)

# User-friendly constructor of composante
composante = function(name) {
    new(Class="Composante", name=name)
}

# Functions get 
# Get the list of models for this composante
setGeneric(
    name="getListModels",
    def=function(object) {standardGeneric("getListModels")}
)

setMethod(
    f="getListModels", 
    signature="Composante",
    definition=function(object) {
        return(object@listModel)
    }
)

# Get the number of models for this composante
setGeneric(
    name="getNbModel",
    def=function(object) {standardGeneric("getNbModel")}
)

setMethod(
    f="getNbModel", 
    signature="Composante",
    definition=function(object) {
        return(object@nbModel)
    }
)

# Get the name of the composante
setGeneric(
    name="getNameComp",
    def=function(object) {standardGeneric("getNameComp")}
)

setMethod(
    f="getNameComp", 
    signature="Composante",
    definition=function(object) {
        return(object@name)
    }
)

# Function to add model(s) in the composante (used by the function in paramList)
setGeneric(
    name="addModel",
    def=function(object, nbToAdd) {standardGeneric("addModel")}
)

setMethod(
    f="addModel",
    signature="Composante",
    definition=function(object, nbToAdd) {
        for(i in 1:nbToAdd) {
            object@nbModel = object@nbModel+1

            newMod = model(object@name, object@nbModel)
            object@listModel = c(object@listModel, newMod)
        }
        rm(newMod)
        validObject(object)
        return(object)
    }
)

# Function to delete model(s) in the composante (used by the function in paramList)
setGeneric(
    name="delModel",
    def=function(object, numModelToDel) {standardGeneric("delModel")}
)

setMethod(
    f="delModel",
    signature="Composante",
    definition=function(object, numModelToDel) {
        compteur = 0
        numModelToDel = sort(numModelToDel)
        # Deletion of the models
        for(i in numModelToDel) {
            object@listModel = object@listModel[-(i-compteur)]
            compteur = compteur+1
        }
        object@nbModel = object@nbModel - length(numModelToDel)
        # Update of the models number
        for(i in 1:object@nbModel) {
            object@listModel[[i]] = setNumModel(object@listModel[[i]], i)
        }
        return(object)
    }
)

# Get the name of the composante
setGeneric(
  name="setComposante",
  def=function(object) {standardGeneric("setComposante")}
)

setMethod(
  f="setComposante", 
  signature="Composante",
  definition=function(object) {
    flag = 0
    while(flag == 0){
      cat("What do you want to do? (press 0 to quit)\n 1: Add a model\n 2: Delete a model\n 3: Change a model")
      scanner = as.numeric(readline())
      if(is.na(scanner) || scanner>3 || scanner<1){
        print("ERROR: Your entry is incorrect, please try again")
      } else if(scanner == 0){
        stop("You stopped the program")
      }else{
        flag = 1
      }
    }
    if(scanner == 1){
      flag = 0
      while(flag == 0){
        cat("How many models do you want to add ?")
        nbToAdd = as.numeric(readline())
        if(is.na(nbToAdd) || nbToAdd<1){
          print("ERROR: Your entry is incorrect, please try again")
        }else{
          flag = 1
        }
      }
      object = addModel(object, nbToAdd)
    } else if(scanner == 2){
      flag = 0
      while(flag == 0){
          if(object@nbModel > 1) {
              print(object)
              cat("Which model do you want to delete ?")
              nbToDel = as.numeric(readline())
              if(is.na(nbToDel) || nbToDel> getNbModel(object)){
                print("ERROR: Your entry is incorrect, please try again")
              }else{
                flag = 1
              }
              object = delModel(object, nbToDel)
          } else {
              stop("There is only one model left in this composante. You can not delete it.")
          }
      }

    } else {
      flag = 0
      while(flag == 0){
        print(object)
        cat("Which model do you want to change ? Please enter the model's number")
        change = as.numeric(readline())
        if(is.na(change) || change> getNbModel(object) || change<1){
          print("ERROR: Your entry is incorrect, please try again")
        }else{
          flag = 1
        }
      }

      

      flag = 0
      while(flag == 0){
        cat("What do you want to change?\n 1. Model function 2. Model parameters")
        choice = as.integer(readline())
        if(is.na(change) || change> getNbModel(object) || choice<1){
          print("ERROR: Your entry is incorrect, please try again")
        }else{
          flag = 1
        }
      }      
      print(object@listModel[[change]])
      if(choice == 1) {
        object@listModel[[change]] = setTypeModel(object@listModel[[change]])
      }
      else if(choice == 2) {
        object@listModel[[change]] = setPrior(object@listModel[[change]])
      }
      
    }
    return(object)
  }
)


# Function to print the parameters of all the composante models
setMethod(
    f="show", 
    signature="Composante",
    definition=function(object) {
        for(i in 1:length(object@listModel)) {
            cat("............... Model:",i,",",object@listModel[[i]]@type_model,"...............\n")
            print(object@listModel[[i]])
        }
    }
)

# Function to assess the prior values for one model
setGeneric(
    name="setResultPriorComp",
    def=function(object, all) {standardGeneric("setResultPriorComp")}
)


setMethod(
    f="setResultPriorComp",
    signature="Composante",
    definition=function(object, all) {
        if(all == 0) { 
            # ask which model the user wants to assess
            cat("[Type 0 to quit] Which model do you want to assess? Type the model number.\n")
            print(object)
            cat("Type", object@nbModel+1, "to assess all models.\n")
            
            flag = -1
            while(flag == -1) {
                choice = as.integer(readline())
                if(choice!=0 && choice<=object@nbModel+1 && choice>0 && !is.na(choice)) {
                    if(choice == object@nbModel+1) {
                        for(i in 1:object@nbModel) {
                            object@listModel[[i]] = setResultPriorMod(object@listModel[[i]],1)
                        }
                    }
                    else {
                        object@listModel[[choice]] = setResultPriorMod(object@listModel[[choice]],0)
                    }
                    flag = 1
                }
                else if (choice == 0 && !is.na(choice)) {
                    stop("Stop the program.")
                }
                else {
                    print("Wrong number, please enter a model number.")
                }
            }
        } else if(all == 1) {

        }
        return(object)
    }
)