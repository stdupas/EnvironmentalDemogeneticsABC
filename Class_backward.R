# Sources for needed classes
source("Class_paramList.R")

# Class backward
setClass(
    Class="Backward", 
    contains="ParamList",
    validity=function(object) {
        # add verification if needed
    }
)

# Constructor of backward
setMethod(
    f="initialize",
    signature="Backward",
    definition=function(.Object) {
        .Object@niche = composante("niche")
        .Object@dispersion = composante("dispersion")
        .Object@mutation = composante("mutation")
        validObject(.Object)
        return(.Object)
    }
)

# User-friendly constructor of backward
backward = function() {
    new(Class="Backward")
}

# Change any thing in the model
setGeneric(
  name="setBackward",
  def=function(object) {standardGeneric("setBackward")}
)

setMethod(
  f="setBackward", 
  signature="Backward",
  definition=function(object) {
    flag = 0
    while(flag == 0){
      cat("Which composante do you want to change? (press 0 to quit)\n 1: Niche\n 2: Dispersion\n 3: Mutation")
      scanner = as.numeric(readline())
      if(is.na(scanner) || scanner>3 || scanner<0){
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
    }
    return(object)
  }
)
