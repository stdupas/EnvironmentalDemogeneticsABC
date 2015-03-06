#Class of paramModel
setClass(
  Class = "ParamModel",
  representation = representation(
    name = "character",
    type_prior = "character",
    param_prior = "numeric",
    result_prior = "list"
    ),
  prototype = prototype(
    name = character(0),
    type_prior = character(0), 
    param_prior = numeric(0),
    result_prior = list(0)
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
    for (i in 1:vec[1]){
      flag = 0
        while(flag == 0){
          print(paste("What do you want for the parameter ", vec[i+1]," ? (press 0 to quit)"))
          scanner = as.numeric(readline())
          if (is.na(scanner) || (i == 1 && scanner<0)){
            print("ERROR: Your entry is incorrect, please try again")
          }else if(scanner == 0){
            stop("You stopped the programm")
          } else {
            flag = 1
          }      
        }
      param = c(param, scanner)
    }
    .Object@param_prior = param
    validObject(.Object)
    return(.Object)
  }
)

#UserFriendly constructor with 1 arguments
paramModel = function(model_num){
  new(Class = "ParamModel", model_num = model_num)
}

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
            #A COMPLETER
            return(object@result_prior)
          }
)


