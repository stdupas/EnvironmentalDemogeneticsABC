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
#Args: object paramList, component name, c(model number, parameter model) or 0 for independant model
setGeneric("getResultPrior",
           function(object,comp,param){standardGeneric("getResultPrior")})

setMethod("getResultPrior", "ParamList",
          function(object,comp,param){
            if(comp == "niche_r") {
                if(param[1] == 0) {
                    return(getResult_prior(object@niche_r@independance))
                } else {
                    return(getResult_prior(object@niche_r@listModel[[param[1]]]@param_model[[param[2]]]))
                }
            } else if(com == "niche_k") {
                if(param[1] == 0) {
                    return(getResult_prior(object@niche_k@independance))
                } else {
                    return(getResult_prior(object@niche_k@listModel[[param[1]]]@param_model[[param[2]]]))
                }
            } else if(comp == "dispersion") {
                return(getResult_prior(object@dispersion@listModel[[param[1]]]@param_model[[param[2]]]))
            } else if(comp == "mutation") {
                return(getResult_prior(object@mutation@listModel[[param[1]]]@param_model[[param[2]]]))
            } else if(comp == "generation") {
                if(param[1] == 0) {
                    return(getResult_prior(object@generation@independance))
                } else {
                    return(getResult_prior(object@generation@listModel[[param[1]]]@param_model[[param[2]]]))
                }
            }
          }
)