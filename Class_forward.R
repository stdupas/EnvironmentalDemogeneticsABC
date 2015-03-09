# Sources for needed classes
source("Class_composante.R")
source("Class_paramList.R")

# Class forward
setClass(
	Class="Forward",
	representation=representation(
		generation="Composante"
	),
	contains="ParamList",
    validity=function(object) {
    	# add verification if needed
    }
)

# Constructor of forward
setMethod(
	f="initialize",
	signature="Forward",
	definition=function(.Object) {
		.Object@niche = composante("niche")
		.Object@dispersion = composante("dispersion")
		.Object@mutation = composante("mutation")
		.Object@generation = composante("generation")
		validObject(.Object)
        return(.Object)
	}
)

# User-friendly constructor of forward
forward = function() {
	new(Class="Forward")
}

# Change any thing in the forward object
setGeneric(
  name="setForward",
  def=function(object) {standardGeneric("setForward")}
)

setMethod(
  f="setForward", 
  signature="Forward",
  definition=function(object) {
    flag = 0
    while(flag == 0){
      cat("Which composante do you want to change? (press 0 to quit)\n 1: Niche\n 2: Dispersion\n 3: Mutation\n 4: Génération")
      scanner = as.numeric(readline())
      if(is.na(scanner) || scanner>4 || scanner<0){
        print("ERROR: Your entry is incorrect, please try again")
      } else if(scanner == 0){
        stop("You stopped the program")
      }else{
        flag = 1
      }
    }
    if(scanner == 1){
      object@niche = setComposante(object@niche)
    } else if (scanner == 2){
      object@dispersion = setComposante(object@dispersion)
    } else if (scanner == 3){
      object@mutation = setComposante(object@mutation)
    } else if (scanner == 4){
      object@generation = setComposante(object@generation)
    }
    return(object)
  }
)