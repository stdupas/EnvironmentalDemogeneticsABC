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
    name="getName",
    def=function(object) {standardGeneric("getName")}
)

setMethod(
    f="getName", 
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

