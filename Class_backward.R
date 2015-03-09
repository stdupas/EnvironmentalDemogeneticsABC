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
      cat("[Type 0 to quit] Which composante do you want to change?\n1: Niche\n2: Dispersion\n3: Mutation")
      scanner = as.numeric(readline())
      if(is.na(scanner) || scanner>3 || scanner<0){
        print("ERROR: Your entry is incorrect, please try again")
      } else if(scanner == 0){
        stop("You have stopped the program")
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
    flag1 = -1
    while(flag1 == -1){
      cat("[Type 0 to quit] What do you want to do? \n1: Assess the prior values of all the components\n2: Assess the prior values of one component\n")
      choice = as.integer(readline())
      if(!is.na(choice) && choice!=0 && choice>0 && choice<=2){
        if(choice == 1){
          object@niche = setResultPriorComp(object@niche, 1)
          object@dispersion = setResultPriorComp(object@dispersion, 1)
          object@mutation = setResultPriorComp(object@mutation, 1)
        } else if (choice == 2){
          flag2 = -1
          while(flag2 == -1) {
            cat("[Type 0 to quit] Which component do you want to assess? \n1: Niche\n2: Dispersion\n3: Mutation\n")
            choice = as.integer(readline())
            if(choice!=0 && choice<=3 && choice>0 && !is.na(choice)) {
              if(choice == 1){
                object@niche = setResultPriorComp(object@niche, 0)
              } else if(choice == 2){
                object@dispersion = setResultPriorComp(object@dispersion, 0)
              } else if(choice == 3){
                object@mutation = setResultPriorComp(object@mutation, 0)
              }
              flag2 = 1
            }
            else if (choice == 0 && !is.na(choice)) {
              stop("You have stopped the program")
            }
            else {
              print("ERROR: Your entry is incorrect, please try again")
            }
          }  
        }
        flag1 = 1
      } else if (choice == 0){
        stop("You have stopped the program")
      } else {
        print("ERROR: Your entry is incorrect, please try again")
      }
    }
    return(object)
  }
)

