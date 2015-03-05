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
    	cat("---------- Composantes : verification ----------\n")
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
        cat("---------- Composantes : initiation ----------\n")
        .Object@name=name
        .Object@nbModel = as.numeric(readline(paste("How many models for",name,"? ")))
        
        mod = NULL
        for(i in 1:.Object@nbModel) {
            mod = c(mod, model(name,i))
        }

        .Object@listModel = mod
        validObject(.Object)
        return(.Object)
    }
)

# User-friendly constructor of composante
composante = function(name) {
    cat("---------- Composantes : construction ----------\n")
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

# Function to add a model in the composante (used by the function in paramList)
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
    }
)
