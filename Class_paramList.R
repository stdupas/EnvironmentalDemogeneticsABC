# Sources for needed classes
source("Class_composante.R")

# Class paramList
setClass(
	Class="ParamList", 
    representation=representation(
    	niche_r="Composante",
        niche_k="Composante",
    	dispersion="Composante",
    	mutation="Composante"
    ),
    validity=function(object) {
    	cat("---------- ParamList : verification ----------\n")
    	# add verification if needed
    }
)

#Function to save the ParamList
setGeneric("saveParamList",
           function(object, file){standardGeneric("saveParamList")})

setMethod("saveParamList", "ParamList",
          function(object, file){
            saveRDS(object, file)
          }
)


#Function to load the ParamList
loadParamList = function(file){
            object = readRDS(file)
            return(object)
}

#Function to get the "result_prior"
setGeneric("getResultPrior",
           function(object,comp,mod,param){standardGeneric("getResultPrior")})

setMethod("getResultPrior", "ParamList",
          function(object,comp,mod,param){
            if(comp == "niche_r") {
                return(getResult_prior(object@niche_r@listModel[[mod]]@param_model[[param]]))
            } else if(com == "niche_k") {
                return(getResult_prior(object@niche_k@listModel[[mod]]@param_model[[param]]))
            } else if(comp == "dispersion") {
                return(getResult_prior(object@dispersion@listModel[[mod]]@param_model[[param]]))
            } else if(comp == "mutation") {
                return(getResult_prior(object@mutation@listModel[[mod]]@param_model[[param]]))
            } else if(comp == "generation") {
                return(getResult_prior(object@generation@listModel[[mod]]@param_model[[param]]))
            }
          }
)