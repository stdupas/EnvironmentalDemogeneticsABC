# Sources for needed classes
source("Class_composante.R")

# Class backward
setClass(
    Class="Backward", 
    contains="ParamList",
    representation=representation(
        method = "character"
    ),
    validity=function(object) {
        # add verification if needed
    }
)

# Constructor of backward
setMethod(
    f="initialize",
    signature="Backward",
    definition=function(.Object, names_stacks) {
        .Object@names_stacks = names_stacks
        flag = 0
        while(flag == 0){
            cat("[Type 0 to quit] What method do you want to use?\n1: ABC\n")#2: Bayesian\n3: Likelihood\n")
            scanner = as.integer(readline())
            #Only for ABC, otherwise change the 1 into a 3
            if(is.na(scanner) || scanner>1 || scanner<0){
                print("ERROR: Your entry is incorrect, please try again")
            } else if(scanner == 0){
                stop("You have stopped the program")
            }else{
                flag = 1
            }
        }
        if (scanner == 1){
            .Object@method = "ABC"
#        } else if (scanner == 2){
 #           .Object@method = "Bayesian"
  #      } else if (scanner == 3){
   #         .Object@method = "Likelihood"
        }
        .Object@niche_r = composante("niche_r", getBackwardMethod(.Object), getBackwardNames_stacks(.Object))
        .Object@niche_k = composante("niche_k", getBackwardMethod(.Object), getBackwardNames_stacks(.Object))
        .Object@dispersion = composante("dispersion", getBackwardMethod(.Object), getBackwardNames_stacks(.Object))
        .Object@mutation = composante("mutation", getBackwardMethod(.Object), getBackwardNames_stacks(.Object))
        validObject(.Object)
        return(.Object)
    }
)

# User-friendly constructor of backward
backward = function(names_stacks) {
    new(Class="Backward", names_stacks = names_stacks)
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
      cat("[Type 0 to quit] Which composante do you want to change?\n1: Niche_r\n2: Niche_k\n3: Dispersion\n4: Mutation")
      scanner = as.integer(readline())
      if(is.na(scanner) || scanner>4 || scanner<0){
        print("ERROR: Your entry is incorrect, please try again")
      } else if(scanner == 0){
        stop("You have stopped the program")
      }else{
        flag = 1
      }
    }
    if(scanner == 1){
      object@niche_r = setComposante(object@niche_r)
    } else if(scanner == 2) {
      object@niche_k = setComposante(object@niche_k)
    } else if (scanner == 3){
      object@dispersion = setComposante(object@dispersion)
    } else if (scanner == 4){
      object@mutation = setComposante(object@mutation)
    }
    return(object)
  }
)

# Get method function
setGeneric(
    name="getBackwardMethod",
    def=function(object) {standardGeneric("getBackwardMethod")}
)

setMethod(
    f = "getBackwardMethod",
    signature = "Backward",
    definition = function(object){
        return(object@method)
    }
)

setGeneric(
    name="getBackwardNames_stacks",
    def=function(object) {standardGeneric("getBackwardNames_stacks")}
)

setMethod(
    f = "getBackwardNames_stacks",
    signature = "Backward",
    definition = function(object){
        return(object@names_stacks)
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
      if(getBackwardMethod(object) == "Bayesian" || getBackwardMethod(object) == "Likelihood"){
          print("ERROR: Your method doesn't allow any prior result")
      } else { 
          # ask which composante the user wants to assess the prior values
          flag1 = -1
          while(flag1 == -1){
              cat("[Type 0 to quit] What do you want to do? \n1: Assess the prior values of all the components\n2: Assess the prior values of one component\n")
              choice = as.integer(readline())
              if(!is.na(choice) && choice!=0 && choice>0 && choice<=2){
                  if(choice == 1){
                      object@niche_r = setResultPriorComp(object@niche_r, 1)
                      object@niche_k = setResultPriorComp(object@niche_k, 1)
                      object@dispersion = setResultPriorComp(object@dispersion, 1)
                      object@mutation = setResultPriorComp(object@mutation, 1)
                  } else if (choice == 2){
                      flag2 = -1
                      while(flag2 == -1) {
                          cat("[Type 0 to quit] Which component do you want to assess? \n1: Niche_r\n2: Niche_k\n3: Dispersion\n4: Mutation\n")
                          choice = as.integer(readline())
                          if(choice!=0 && choice<=3 && choice>0 && !is.na(choice)) {
                              if(choice == 1){
                                  object@niche_r = setResultPriorComp(object@niche_r, 0)
                              } else if(choice == 2) {
                                  object@niche_k = setResultPriorComp(object@niche_k, 0)
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
      }
    return(object)
  }
)

