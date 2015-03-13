############### CLASS: "ParamModel" ###############

############### Creation of class "ParamModel" ###############
    # method: a string containing the method used by the paramList (ABC, Bayesian, Likelihood)
    # name: a vector containing ("model:", number of the model)
    # type_prior: a string containing the function used for the prior
    # param_prio: a list containing the value of each parameter
    # param_name: a vector containing the names of the parameters
    # result_prior: a vector containing the priors values
setClass(
    Class = "ParamModel",
    representation = representation(
        method = "character",
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
            stop("[ ParamModel : verification ] no prior type given")
        } else{      
            return (TRUE)
        }
    }
    
)


############### Initiator ###############
    # Function to initialize an object of class "ParamModel"
    # Args:
    #       model_num: an integer containing the number of the actual model
    #       method: a string containing the method used by the paramList (ABC, Bayesian, Likelihood)
setMethod(
    f = "initialize",
    signature = "ParamModel",
    definition = function(.Object, model_num, method){
        .Object@name=c("model:",model_num)
        .Object@method = method
        print("[Type 0 to quit] What function do you want to use for prior ?")
        data_fct = read.table("functions.txt", sep = ";", header = TRUE, as.is=rep(TRUE, 4))
        if(method == "Likelihood"){
            possible = which(data_fct[,1]=="prior" & data_fct[,2] == "uniform")
        } else {
            possible = which(data_fct[,1]=="prior")
        }
        num = c(1:length(possible))
        aff = rep(":", length(possible))
        flag =0
        while(flag ==0){
            print(paste(num, aff,data_fct[possible,2]))    
            scanner = as.integer(readline())
            if (is.na(scanner) || (scanner<0 || scanner > length(possible))){
                print("ERROR: Your entry is incorrect, please try again")
            }else if(scanner == 0){
                stop("You have stopped the programm")
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
                if (i == 1){
                    scanner = as.integer(readline())
                } else {
                    scanner = as.numeric(readline())
                }
                if (is.na(scanner) || (i == 1 && scanner<1)){
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

    # User-friendly constructor
paramModel = function(model_num, method ){
    new(Class = "ParamModel", model_num = model_num, method = method)
}


#############################################
############### Get functions ###############
#############################################

############### getType_Prior ###############
    # Function to get the "type_prior" attribut (prior function)
setGeneric("getType_prior",
           function(object){standardGeneric("getType_prior")})

setMethod("getType_prior", "ParamModel",
          function(object){
              return(object@type_prior)
          }
)

############### getParam_Prior ###############
    # Function to get the "param_prior" attribut
setGeneric("getParam_prior",
           function(object){standardGeneric("getParam_prior")})

setMethod("getParam_prior", "ParamModel",
          function(object){
              return(object@param_prior)
          }
)

############### getNameParamModel ###############
    # Function to get the "name" attribut
setGeneric("getNameParamModel",
           function(object){standardGeneric("getNameParamModel")})

setMethod("getNameParamModel", "ParamModel",
          function(object){
              return(object@name)
          }
)

############### getResult_Prior ###############
    # Function to get the "result_prior" attribut
setGeneric("getResult_prior",
           function(object){standardGeneric("getResult_prior")})

setMethod("getResult_prior", "ParamModel",
          function(object){
              return(object@result_prior)
          }
)

############### getParamMethod ###############
    # Function to get the "method" attribut
setGeneric("getParamMethod",
           function(object){standardGeneric("getParamMethod")})

setMethod("getParamMethod", "ParamModel",
          function(object){
              return(object@method)
          }
)

############### getParam_name ###############
    # Function to get the "param_name" attribut
setGeneric("getParam_name",
           function(object){standardGeneric("getParam_name")})

setMethod("getParam_name", "ParamModel",
          function(object){
              return(object@param_name)
          }
)


#############################################
############### Set functions ###############
#############################################

############### setType_prior ###############
    # Function to update type_prior
setGeneric(
    name="setType_prior",
    def=function(object) {standardGeneric("setType_prior")}
)

setMethod(
    f="setType_prior",
    signature="ParamModel",
    definition=function(object) {
        if(getParamMethod(object) == "Likelihood"){
            print("You can not change the prior type with the method: Likelihood")
            return(object)
        } else {
            print(paste("The actual prior function is: ", getType_prior(object)))
            newObject = new(Class = "ParamModel", model_num = getNameParamModel(object)[2], method = getParamMethod(object))
            return(newObject)   
        }
    }
)

############### setParam_prior ###############
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
        cat(paste(num,":",getParam_name(object),"=",getParam_prior(object),"\n"))
        flag = 0
        while(flag == 0){
            cat("[Type 0 to quit] Which one do you want to change?")
            scanner = as.integer(readline())
            if (is.na(scanner) || scanner>length(getParam_prior(object)) || scanner<0){
                print("ERROR: Your entry is incorrect, please try again")
            } else if(scanner == 0){
                stop("You have stopped the program")
            }else{
                flag = 1
            }      
        }
        flag = 0
        while(flag == 0){
            cat("What is its new valor?")
            if(scanner ==1){
                new_val = as.integer(readline())
            } else {
                new_val = as.numeric(readline())
            }
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

############### setResult_prior ###############
    # Function to change the Result_prior
setGeneric("setResult_prior",
           function(object){standardGeneric("setResult_prior")})

setMethod("setResult_prior", "ParamModel",
          function(object){
              object@result_prior = do.call(getType_prior(object), as.list(getParam_prior(object)))
              return(object)
          }
)

############### delResult_prior ###############
    # Function to clean the Result_prior
setGeneric("delResult_prior",
           function(object){standardGeneric("delResult_prior")})

setMethod("delResult_prior", "ParamModel",
          function(object){
              object@result_prior = numeric(0)
              return(object)
          }
)


#############################################
############### Show functions ##############
#############################################

############### show ###############
    # Function to print the parameters of the prior functions
setMethod(
    f="show", 
    signature="ParamModel",
    definition=function(object) {
        param = paste(getParam_name(object),"=",getParam_prior(object))
        cat(getType_prior(object),"(",param,")\n")
    }
)


#############################################
############## Static functions #############
#############################################

############### findFunctionFromFile ###############
    # Function to find a given function in the file "functions.txt"
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

