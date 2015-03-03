# Sources for needed classes
source("ClassModel.R")

# Class composante
setClass(
	Class="composante", 
    representation=representation(
    	listModel="Model",
    	nbModel="numeric"
    ),
    prototype=prototype(
    	listModel=NULL,
    	nbModel=1
    ),
    validity=function(object) {
    	cat("---------- Composantes : verification ----------")
    	if(object@nbModel <= 0) {
    		stop("[Composante validation] Number of models is not positive")
    	}
    	if(length(object@listModel)!=object@nbModel) {
    		stop("[Composante validation] Number of models is different from nbModel")
    	}
    	return(TRUE)
    }
)