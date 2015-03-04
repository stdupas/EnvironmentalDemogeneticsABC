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
            mod = c(mod, model())
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
