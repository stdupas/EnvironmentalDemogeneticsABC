# Sources for needed classes
source("Class_composante.R")

# Class forward
setClass(
    Class="Forward",
    representation=representation(
        method = "character",
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
    definition=function(.Object, names_stacks) {
        .Object@names_stacks = names_stacks
        flag = 0
        while(flag == 0){
            cat("[Type 0 to quit] What method do you want to use?\n1: Bayesian\n2: Likelihood\n")
            scanner = as.integer(readline())
            if(is.na(scanner) || scanner>2 || scanner<0){
                print("ERROR: Your entry is incorrect, please try again")
            } else if(scanner == 0){
                stop("You have stopped the program")
            }else{
                flag = 1
            }
        }
        if (scanner == 1){
            .Object@method = "Bayesian"
        } else if (scanner == 2){
            .Object@method = "Likelihood"
        } 
        .Object@niche_r = composante("niche_r", getForwardMethod(.Object), getForwardNames_stacks(.Object))
        .Object@niche_k = composante("niche_k", getForwardMethod(.Object), getForwardNames_stacks(.Object))
        .Object@dispersion = composante("dispersion", getForwardMethod(.Object), getForwardNames_stacks(.Object))
        .Object@mutation = composante("mutation", getForwardMethod(.Object), getForwardNames_stacks(.Object))
        .Object@generation = composante("generation", getForwardMethod(.Object), getForwardNames_stacks(.Object))
        validObject(.Object)
        return(.Object)
    }
)

# User-friendly constructor of forward
forward = function(names_stacks) {
    new(Class="Forward", names_stacks = names_stacks)
}

# Change any thing in the forward object
setGeneric(
  name="setForward",
  def=function(object) {standardGeneric("setForward")}
)

# Get method function

setMethod(
  f="setForward", 
  signature="Forward",
  definition=function(object) {
    flag = 0
    while(flag == 0){
      cat("[Type 0 to quit] Which composante do you want to change?\n1: Niche_r\n2: Niche_k\n3: Dispersion\n4: Mutation\n5: Génération")
      scanner = as.integer(readline())
      if(is.na(scanner) || scanner>5 || scanner<0){
        print("ERROR: Your entry is incorrect, please try again")
      } else if(scanner == 0){
        stop("You stopped the program")
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
    } else if (scanner == 5){
      object@generation = setComposante(object@generation)
    }
    return(object)
  }
)

# Get method function
setGeneric(
    name="getForwardMethod",
    def=function(object) {standardGeneric("getForwardMethod")}
)

setMethod(
    f = "getForwardMethod",
    signature = "Forward",
    definition = function(object){
        return(object@method)
    }
)


setGeneric(
    name="getForwardNames_stacks",
    def=function(object) {standardGeneric("getForwardNames_stacks")}
)

setMethod(
    f = "getForwardNames_stacks",
    signature = "Forward",
    definition = function(object){
        return(object@names_stacks)
    }
)


# Function to assess the prior values
setGeneric(
  name="setResultPriorFor",
  def=function(object) {standardGeneric("setResultPriorFor")}
)


setMethod(
  f="setResultPriorFor",
  signature="Forward",
  definition=function(object) {
      if(getForwardMethod(object) == "Bayesian" || getForwardMethod(object) == "Likelihood"){
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
                      object@generation = setResultPriorComp(object@generation, 1)
                  } else if (choice == 2){
                      flag2 = -1
                      while(flag2 == -1) {
                          cat("[Type 0 to quit] Which component do you want to assess? \n1: Niche_r\n2: Niche_k\n3: Dispersion\n4: Mutation\n5: Generation\n")
                          choice = as.integer(readline())
                          if(choice!=0 && choice<=4&& choice>0 && !is.na(choice)) {
                              if(choice == 1){
                                  object@niche_r = setResultPriorComp(object@niche_r, 0)
                              } else if(choice == 2) {
                                  object@niche_k = setResultPriorComp(object@niche_k, 0)
                              } else if(choice == 2){
                                  object@dispersion = setResultPriorComp(object@dispersion, 0)
                              } else if(choice == 3){
                                  object@mutation = setResultPriorComp(object@mutation, 0)
                              } else if(choice == 4){
                                  object@generation = setResultPriorComp(object@generation, 0)
                              }
                              flag2 = 1
                          }
                          else if (choice == 0 && !is.na(choice)) {
                              stop("You have stopped the program")
                          } else {
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