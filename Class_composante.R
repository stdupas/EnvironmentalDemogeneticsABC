# Sources for needed classes
source("Class_model.R")

# Class composante
setClass(
    Class="Composante", 
    representation=representation(
        method = "character",
        name="character",
        listModel="list",
        nbModel="numeric",
        type_combinaison = "character",
        independance = "ParamModel"
    ),
    validity=function(object) {
        if(object@nbModel <= 0) {
            stop("[Composante validation] Number of models is not positive\n")
        }
        if(length(object@listModel)!=object@nbModel) {
            stop("[Composante validation] Number of models is different from nbModel\n")
        }
        return(TRUE)
    }
)


setMethod(
    f="initialize",
    signature="Composante",
    definition=function(.Object, name, method) {
        .Object@name=name
        .Object@method=method
        flag = -1
        while(flag == -1){
            cat("What is the combinaison method for the component ", .Object@name," ?\n1: Additive\n2: Multiplicative\n")
            scanner = as.integer(readline())
            if(!is.na(scanner) && scanner>0 && scanner<3) {
                if(scanner == 1){
                    .Object@type_combinaison = "Additive"
                } else if(scanner == 2){
                    .Object@type_combinaison = "Multiplicative"
                }
                flag = 1
            } else {
                print("ERROR: Your entry is incorrect, please try again")
            } 
        }
        if(.Object@name == "niche_r" || .Object@name == "niche_k" || .Object@name == "generation"){
            cat("=========== Creation of the independant valor ==========\n")
            .Object@independance = paramModel(0, .Object@method)
        }

        # repeat while the given number is incorrect
        flag = -1
        while(flag == -1) {          
            choice_number = as.integer(readline(paste("[Type 0 to quit] How many models for the component",name,"?")))
            if(choice_number!=0 && choice_number>0 && !is.na(choice_number)) {
                .Object@nbModel = choice_number
                flag = 1
            }
            else if(choice_number==0 && !is.na(choice_number)) {
                stop("You have stopped the program")
            }
            else {
                print("ERROR: Your entry is incorrect, please try again")
            } 
        }
        
        mod = NULL
        for(i in 1:.Object@nbModel) {
            print(paste("========== Composante : ",name,", model n°", i," =========="))
            mod = c(mod, model(name,i, .Object@method))
        }
        
        .Object@listModel = mod
        validObject(.Object)
        return(.Object)
    }
)

# User-friendly constructor of composante
composante = function(name, method) {
    new(Class="Composante", name=name, method = method)
}

# Functions get 
# Get the list of models for this composante
setGeneric(
    name="getListModels",
    def=function(object) {standardGeneric("getListModels")}
)

setMethod(
    f="getListModels", 
    signature="Composante",
    definition=function(object) {
        return(object@listModel)
    }
)

# Get the number of models for this composante
setGeneric(
    name="getNbModel",
    def=function(object) {standardGeneric("getNbModel")}
)

setMethod(
    f="getNbModel", 
    signature="Composante",
    definition=function(object) {
        return(object@nbModel)
    }
)

# Get the name of the composante
setGeneric(
    name="getNameComp",
    def=function(object) {standardGeneric("getNameComp")}
)

setMethod(
    f="getNameComp", 
    signature="Composante",
    definition=function(object) {
        return(object@name)
    }
)

# Function to add model(s) in the composante (used by the function in paramList)
setGeneric(
    name="addModel",
    def=function(object, nbToAdd) {standardGeneric("addModel")}
)

setMethod(
    f="addModel",
    signature="Composante",
    definition=function(object, nbToAdd) {
        for(i in 1:nbToAdd) {
            object@nbModel = object@nbModel+1
            
            newMod = model(object@name, object@nbModel)
            object@listModel = c(object@listModel, newMod)
        }
        rm(newMod)
        validObject(object)
        return(object)
    }
)

# Function to delete model(s) in the composante (used by the function in paramList)
setGeneric(
    name="delModel",
    def=function(object, numModelToDel) {standardGeneric("delModel")}
)

setMethod(
    f="delModel",
    signature="Composante",
    definition=function(object, numModelToDel) {
        compteur = 0
        numModelToDel = sort(numModelToDel)
        # Deletion of the models
        for(i in numModelToDel) {
            object@listModel = object@listModel[-(i-compteur)]
            compteur = compteur+1
        }
        object@nbModel = object@nbModel - length(numModelToDel)
        # Update of the models number
        for(i in 1:object@nbModel) {
            object@listModel[[i]] = setNumModel(object@listModel[[i]], i)
        }
        return(object)
    }
)

# Function that allows the user to ad, delete or modify a model
setGeneric(
    name="setComposante",
    def=function(object) {standardGeneric("setComposante")}
)

setMethod(
    f="setComposante", 
    signature="Composante",
    definition=function(object) {
        flag = 0
        while(flag == 0){
            if(object@name == "niche_k" || object@name == "niche_r" || object@name == "generation") {
                cat("[Type 0 to quit] What do you want to do?\n1: Add a model\n2: Delete a model\n3: Change a model\n4: Change the independant model")
            } else {
                cat("[Type 0 to quit] What do you want to do?\n1: Add a model\n2: Delete a model\n3: Change a model")
            }
            scanner = as.numeric(readline())
            if(is.na(scanner) || scanner>5 || scanner<0){
                print("ERROR: Your entry is incorrect, please try again")
            } else if(scanner == 0){
                stop("You have stopped the program")
            }else{
                if(scanner == 4 && object@name != "niche_k" && object@name != "niche_r" && object@name != "generation") {
                    print("ERROR: Your entry is incorrect, please try again")
                } else {
                    flag = 1
                }
            }
        }
        if(scanner == 1){
            flag = 0
            while(flag == 0){
                cat("How many models do you want to add ?")
                nbToAdd = as.numeric(readline())
                if(is.na(nbToAdd) || nbToAdd<1){
                    print("ERROR: Your entry is incorrect, please try again")
                }else{
                    flag = 1
                    object = addModel(object, nbToAdd)
                }
            }
        } else if(scanner == 2){
            flag = 0
            while(flag == 0){
                if(object@nbModel > 1) {
                    print(object)
                    cat("Which model do you want to delete ?")
                    nbToDel = as.numeric(readline())
                    if(is.na(nbToDel) || nbToDel> getNbModel(object)){
                        print("ERROR: Your entry is incorrect, please try again")
                    }else{
                        flag = 1
                        object = delModel(object, nbToDel)
                    }
                } else {
                    stop("There is only one model left in this composante. You can not delete it.")
                }
            }
            
        } else if(scanner == 3) {
            flag = 0
            while(flag == 0){
                print(object)
                cat("Which model do you want to change ? Please enter the model's number")
                change = as.numeric(readline())
                if(is.na(change) || change> getNbModel(object) || change<1){
                    print("ERROR: Your entry is incorrect, please try again")
                }else{
                    flag = 1
                }
            }
            flag = 0
            while(flag == 0){
                cat("What do you want to change?\n 1. Model function 2. Model parameters")
                choice = as.integer(readline())
                if(is.na(change) || change> getNbModel(object) || choice<1){
                    print("ERROR: Your entry is incorrect, please try again")
                }else{
                    flag = 1
                }
            }      
            print(object@listModel[[change]])
            if(choice == 1) {
                object@listModel[[change]] = setTypeModel(object@listModel[[change]])
            }
            else if(choice == 2) {
                object@listModel[[change]] = setPrior(object@listModel[[change]])
            }
        } else if (scanner == 4){
            flag2 = -1
            while(flag2 == -1) {
                # Ask if the user wants to change the prior function or only a parameter of a prior function
                cat("[Type 0 to quit] Change the prior function or the value of a parameter value ?\n")
                cat("1: Prior function\n2: Parameter value")
                choice_number2 = as.integer(readline())
                if(choice_number2!=0 && choice_number2<=2 && choice_number2>0 && !is.na(choice_number2)) {
                    if(choice_number2==1) {
                        object@independance = setType_prior(object@independance)
                    }
                    else if(choice_number2==2) {
                        object@independance = setParam_prior(object@independance)
                    }                    
                    flag2 = 1
                }
                else if(choice_number2==0 && !is.na(choice_number2)) {
                    stop("You have stopped the program")
                }else {
                    print("ERROR: Your entry is incorrect, please try again")
                }
            }
            flag = 1
        } else if(choice_number==0 && !is.na(choice_number)) {
            stop("You have stopped the program")
        }
        else {
            print("ERROR: Your entry is incorrect, please try again")
        }
        return(object)
    }
)


# Function to print the parameters of all the composante models
setMethod(
    f="show", 
    signature="Composante",
    definition=function(object) {
        cat(".................... Combination type ....................\n")
        cat("           ")
        cat(object@type_combinaison,"\n")
        if(object@name == "niche_r" || object@name == "niche_k" || object@name == "generation") {
            cat(".................... Independant model ....................\n")
            cat("           ")
            print(object@independance)
        }
        for(i in 1:length(object@listModel)) {
            cat("............... Model:",i,",",object@listModel[[i]]@type_model,"...............\n")
            print(object@listModel[[i]])
        }
    }
)

# Function to assess the prior values for one model
setGeneric(
    name="setResultPriorComp",
    def=function(object, all) {standardGeneric("setResultPriorComp")}
)


setMethod(
    f="setResultPriorComp",
    signature="Composante",
    definition=function(object, all) {
        # if there is an independant model
        if(object@name == "niche_k" || object@name == "niche_r" || object@name == "generation") {
            if(all == 0) { 
                # ask which model the user wants to assess
                cat("[Type 0 to quit] Which model do you want to assess? Type the model number.\n")
                print(object)
                cat("Type", object@nbModel+1, "to assess the independant model.\n")
                cat("Type", object@nbModel+2, "to assess all models.\n")
                

                flag = -1
                while(flag == -1) {
                    choice = as.integer(readline())
                    if(choice!=0 && choice<=object@nbModel+2 && choice>0 && !is.na(choice)) {
                        if(choice == object@nbModel+1) {
                            object@independance = setResult_prior(object@independance)
                        } else if(choice == object@nbModel+2) {
                            object@independance = setResult_prior(object@independance)
                            for(i in 1:object@nbModel) {
                                object@listModel[[i]] = setResultPriorMod(object@listModel[[i]],1)
                            }
                        }
                        else {
                            object@listModel[[choice]] = setResultPriorMod(object@listModel[[choice]],0)
                        }
                        flag = 1
                    }
                    else if (choice == 0 && !is.na(choice)) {
                        stop("You have stopped the program")
                    }
                    else {
                        print("ERROR: Your entry is incorrect, please try again")
                    }
                }
            } else if(all == 1) {
                object@independance = setResult_prior(object@independance)
                for(i in 1:object@nbModel) {
                    object@listModel[[i]] = setResultPriorMod(object@listModel[[i]],1)
                }
            }

        # if there is no independant model            
        } else {
            if(all == 0) { 
                # ask which model the user wants to assess
                cat("[Type 0 to quit] Which model do you want to assess? Type the model number.\n")
                print(object)
                cat("Type", object@nbModel+1, "to assess all models.\n")
                

                flag = -1
                while(flag == -1) {
                    choice = as.integer(readline())
                    if(choice!=0 && choice<=object@nbModel+1 && choice>0 && !is.na(choice)) {
                        if(choice == object@nbModel+1) {
                            for(i in 1:object@nbModel) {
                                object@listModel[[i]] = setResultPriorMod(object@listModel[[i]],1)
                            }
                        }
                        else {
                            object@listModel[[choice]] = setResultPriorMod(object@listModel[[choice]],0)
                        }
                        flag = 1
                    }
                    else if (choice == 0 && !is.na(choice)) {
                        stop("You have stopped the program")
                    }
                    else {
                        print("ERROR: Your entry is incorrect, please try again")
                    }
                }
            } else if(all == 1) {
                for(i in 1:object@nbModel) {
                    object@listModel[[i]] = setResultPriorMod(object@listModel[[i]],1)
                }
            }
        }
        return(object)
    }
)