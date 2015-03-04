# Sources for needed classes
source("ClassModel.R")

# Class composante
setClass(
	Class="Composante", 
    representation=representation(
    	listModel="list",
    	nbModel="numeric"
    ),
    prototype=prototype(
    	listModel=NULL,
    	nbModel=1
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
setMethod(
    f="initialize",
    signature="Composante",
    definition=function(.Object, listModel, nbModel) {
        cat("---------- Composantes : initiation ----------\n")
        if(!missing(listModel) && !missing(nbModel)) {
            .Object@listModel = listModel
            .Object@nbModel = nbModel
            validObject(.Object)
        }
        else {
            stop("[Composante initiation] Missing argument(s)\n")
        }
        return(.Object)
    }
)

# User-friendly constructor of composante
composante = function(listModel, nbModel) {
    cat("---------- Composantes : construction ----------\n")
    new(Class="Composante", listModel=listModel, nbModel=nbModel)
}
