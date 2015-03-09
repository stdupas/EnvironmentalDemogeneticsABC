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

# Function to assess the prior values
setGeneric(
  name="setResultPriorBack",
  def=function(object) {standardGeneric("setResultPriorBack")}
)


setMethod(
  f="setResultPriorBack",
  signature="Backward",
  definition=function(object) {
    # ask which composante the user wants to assess the prior values
    
    flag = -1
    while(flag == -1) {
      cat("[Type 0 to quit] Which component do you want to assess? \n1: Niche\n2: Dispersion\n3: Mutation\n")
      choice = as.integer(readline())
      if(choice!=0 && choice<=3 && choice>0 && !is.na(choice)) {
        if(choice == 1){
          object@niche = setResult_priorComp(object@niche)
          flag = 1
        } else if(choice == 2){
          object@dispersion = setResult_priorComp(object@dispersion)
          flag = 1
        } else if(choice == 3){
          object@generation = setResult_priorComp(object@genertion)
          flag = 1
        }
      }
      else if (choice == 0 && !is.na(choice)) {
        stop("Stop the program.")
      }
      else {
        print("Wrong number, please type a number in the list :")
      }
    }
    return(object)
  }
)

