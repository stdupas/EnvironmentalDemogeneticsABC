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
        if(object@method == "Bayesian" || object@method == "Likelihood") {
            cat("ERROR: Your method doesn't allow any prior result")
        } else if(comp == "niche_r") {
            if(param[1] == 0) {
                res = getResult_prior(object@niche_r@independance)
            } else {
                res = getResult_prior(object@niche_r@listModel[[param[1]]]@param_model[[param[2]]])
            }
        } else if(comp == "niche_k") {
            if(param[1] == 0) {
                res = getResult_prior(object@niche_k@independance)
            } else {
                res = getResult_prior(object@niche_k@listModel[[param[1]]]@param_model[[param[2]]])
            }
        } else if(comp == "dispersion") {
            res = getResult_prior(object@dispersion@listModel[[param[1]]]@param_model[[param[2]]])
        } else if(comp == "mutation") {
            res = getResult_prior(object@mutation@listModel[[param[1]]]@param_model[[param[2]]])
        } else if(comp == "generation") {
            if(param[1] == 0) {
                res = getResult_prior(object@generation@independance)
            } else {
                res = getResult_prior(object@generation@listModel[[param[1]]]@param_model[[param[2]]])
            }
        }
        if(length(res) == 0) {
            cat("ERROR: You did not compute the result prior yet.")
        } else {
            return(res)
        }
    }
)