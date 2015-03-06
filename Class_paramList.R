# Sources for needed classes
source("Class_composante.R")

# Class paramList
setClass(
	Class="ParamList", 
    representation=representation(
    	niche="Composante",
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

