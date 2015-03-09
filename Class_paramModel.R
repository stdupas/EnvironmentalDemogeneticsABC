#Class of paramModel
setClass(
  Class = "ParamModel",
  representation = representation(
    name = "character",
    type_prior = "character",
    param_prior = "numeric",
    param_name = "character",
    result_prior = "numeric"
    ),
  prototype = prototype(
    name = character(0),
    type_prior = character(0), 
    param_prior = numeric(0),
    param_name = character(0),
    result_prior = numeric(0)
    ),
  validity = function(object){
    if((is.null(object@type_prior))){
      stop("[ ParamModel : verificatiteston ] no prior type given")
    } else{      
      return (TRUE)
    }
  }
  
)


#Initiateur with 1 arguments
setMethod(
  f = "initialize",
  signature = "ParamModel",
  definition = function(.Object, model_num){
    .Object@name=c("model:",model_num)
    print("What function do you want to use for prior ? (press 0 to quit)")
    data_fct = read.table("functions.txt", sep = ";", header = TRUE, as.is=rep(TRUE, 4))
    possible = which(data_fct[,1]=="prior")
    num = c(1:length(possible))
    aff = rep(":", length(possible))
    flag =0
    while(flag ==0){
      print(paste(num, aff,data_fct[possible,2]))    
      scanner = as.numeric(readline())
      if (is.na(scanner) || (scanner<0 || scanner > length(possible))){
        print("ERROR: Your entry is incorrect, please try again")
      }else if(scanner == 0){
        stop("You stopped the programm")
      } else {
        flag = 1
      }      
    }
    .Object@type_prior = data_fct[possible[scanner],2]
    vec = findFunctionFromFile("prior", .Object@type_prior)
    param=NULL
    param_name = NULL
    for (i in 1:vec[1]){
      flag = 0
        while(flag == 0){
          print(paste("What do you want for the hyper-parameter ", vec[i+1]," ?"))
          scanner = as.numeric(readline())
          if (is.na(scanner) || (i == 1 && scanner<0)){
            print("ERROR: Your entry is incorrect, please try again")
          } else {
            flag = 1
          }      
        }
      param = c(param, scanner)
      param_name = c(param_name, vec[i+1])
    }
    .Object@param_prior = param
    .Object@param_name = param_name
    validObject(.Object)
    return(.Object)
  }
)

#UserFriendly constructor with 1 arguments
paramModel = function(model_num){
  new(Class = "ParamModel", model_num = model_num)
}



######################### METHODS ##################################

#Function to get the "type_prior" attribut
setGeneric("getType_prior",
           function(object){standardGeneric("getType_prior")})

setMethod("getType_prior", "ParamModel",
          function(object){
            return(object@type_prior)
          }
)

#Function to get the "param_prior" attribut
setGeneric("getParam_prior",
           function(object){standardGeneric("getParam_prior")})

setMethod("getParam_prior", "ParamModel",
          function(object){
            return(object@param_prior)
          }
)

#Function to get the "name" attribut
setGeneric("getNameParamModel",
           function(object){standardGeneric("getNameParamModel")})

setMethod("getNameParamModel", "ParamModel",
          function(object){
            return(object@name)
          }
)

#Function to get the "param_prior" attribut
setGeneric("getParam_prior",
           function(object){standardGeneric("getParam_prior")})

setMethod("getParam_prior", "ParamModel",
          function(object){
            return(object@param_prior)
          }
)

#Function to get the "result_prior" attribut
setGeneric("getResult_prior",
           function(object){standardGeneric("getResult_prior")})

setMethod("getResult_prior", "ParamModel",
          function(object){
            return(object@result_prior)
          }
)


# Function to update type_prior
setGeneric(
  name="setType_prior",
  def=function(object) {standardGeneric("setType_prior")}
)

setMethod(
  f="setType_prior",
  signature="ParamModel",
  definition=function(object) {
    print(paste("The actual prior function is: ", getType_prior(object)))
    newObject = new(Class = "ParamModel", model_num = getNameParamModel(object)[2])
    return(newObject)
  }
)

# Function to update param_prior
setGeneric(
  name="setParam_prior",
  def=function(object) {standardGeneric("setParam_prior")}
)

setMethod(
  f="setParam_prior",
  signature="ParamModel",
  definition=function(object) {
    print("The actual prior hyper-parameters are: ")
    num = c(1:length(getParam_prior(object)))
    cat(paste(num,":",object@param_name,"\n"))
    flag = 0
    while(flag == 0){
      cat("Which one do you want to change? (press 0 to quit)")
      scanner = as.numeric(readline())
      if (is.na(scanner) || (scanner>length(getParam_prior(object)) && scanner<0)){
        print("ERROR: Your entry is incorrect, please try again")
      } else if(scanner == 0){
        stop("You stopped the program")
      }else{
        flag = 1
      }      
    }
    flag = 0
    while(flag == 0){
      cat("What is its new valor?")
      new_val = as.numeric(readline())
      if (is.na(new_val) || (scanner == 1 && new_val<=0)){
        print("ERROR: Your entry is incorrect, please try again")
      } else {
        flag = 1
      }      
    }
    object@param_prior[scanner] = new_val  
    return(object)
  }
)



findFunctionFromFile = function(model_type,fct_name){
  data_fct = read.table("functions.txt", sep = ";", header = TRUE)
  if(!is.element(fct_name, data_fct[,2])){
    stop("The function: ", fct_name," does not match any knowm function")
  } else {
    indice = match(fct_name, data_fct[,2])
    if(data_fct[indice,1]!=model_type){
      stop("This function is not a ", model_type, " function")
    } else {
      paul=unlist(strsplit((toString(data_fct[indice,4])),","))
      return (c(data_fct[indice,3], as.vector(paul)))
    }
  }
}

#Function to change the Result_prior
setGeneric("setResult_prior",
           function(object){standardGeneric("setResult_prior")})

setMethod("setResult_prior", "ParamModel",
          function(object){
            object@result_prior = do.call(getType_prior(object), as.list(getParam_prior(object)))
            return(object)
          }
)

#Function to clean the Result_prior
setGeneric("delResult_prior",
           function(object){standardGeneric("delResult_prior")})

setMethod("delResult_prior", "ParamModel",
          function(object){
            object@result_prior = numeric(0)
            return(object)
          }
)


# Function to print the parameters of the prior functions
setMethod(
    f="show", 
    signature="ParamModel",
    definition=function(object) {
        param = paste(object@param_name,"=",object@param_prior)
        cat(object@type_prior,"(",param,")\n")
    }
)
