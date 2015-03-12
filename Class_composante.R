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
        independance = "ParamModel",
        nb_stacks = "numeric",
        names_stacks = "character",
        remain_layers = "character"
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
    definition=function(.Object, name, method, nb_stacks, names_stacks) {
        .Object@name=name
        .Object@method=method
        .Object@nb_stacks = nb_stacks
        .Object@names_stacks = names_stacks
        flag = -1
        while(flag == -1){
            cat("What is the combinaison method for the component ", name," ?\n1: Additive\n2: Multiplicative\n")
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
        if(name == "niche_r" || name == "niche_k" || name == "generation"){
            cat("=========== Creation of the independant valor ==========\n")
            .Object@independance = paramModel(0, getMethodComp(.Object))
        }
        
        if(name == "niche_r" || name == "niche_k" || name == "generation"){
            vec2 = getNames_stacks(.Object)
            vec3 = NULL
            compteur = 0
            sortie = 0
            mod = NULL
            while (compteur < getNb_stacks(.Object) && sortie == 0){
                vec = 1:(getNb_stacks(.Object) - length(vec3))
                cat("Which variable do you wan't to build ? (type 0 if you don't want to build any other model)\n")
                cat(paste(vec,":", vec2),"\n")
                scanner = as.integer(readline())
                if((scanner<=length(vec) && scanner>0) && !is.na(scanner)) {
                    vec3 = c(vec3, scanner)
                    .Object@nbModel = length(vec3)
                    print(.Object@nbModel)
                    compteur = compteur + 1                    
                    mod = c(mod, model(name,compteur, getMethodComp(.Object), vec2[vec3[compteur]]))
                    vec2 = vec2[-vec3[compteur]]
                }
                else if(scanner==0 && !is.na(scanner)) {
                    sortie = 1
                }
                else {
                    print("ERROR: Your entry is incorrect, please try again")
                } 
            }
            .Object@remain_layers = vec2
        } else {
            
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
            for(i in 1:getNbModel(.Object)) {
                print(paste("========== Composante : ",name,", model n°", i," =========="))
                mod = c(mod, model(name,i, getMethodComp(.Object), "STANDARD"))
            }
        }
        
        .Object@listModel = mod
        validObject(.Object)
        return(.Object)
    }
)


# User-friendly constructor of composante
composante = function(name, method, nb_stacks, names_stacks) {
    new(Class="Composante", name=name, method = method, nb_stacks = nb_stacks, names_stacks = names_stacks)
}

#################################### GET METHODS ####################################

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

# Get the combinaison type of the composante
setGeneric(
    name="getType_combinaison",
    def=function(object) {standardGeneric("getType_combinaison")}
)

setMethod(
    f="getType_combinaison", 
    signature="Composante",
    definition=function(object) {
        return(object@type_combinaison)
    }
)

# Get the method of the composante
setGeneric(
    name="getMethodComp",
    def=function(object) {standardGeneric("getMethodComp")}
)

setMethod(
    f="getMethodComp", 
    signature="Composante",
    definition=function(object) {
        return(object@method)
    }
)

setGeneric(
    name="getNb_stacks",
    def=function(object) {standardGeneric("getNb_stacks")}
)

setMethod(
    f="getNb_stacks", 
    signature="Composante",
    definition=function(object) {
        return(object@nb_stacks)
    }
)

setGeneric(
    name="getNames_stacks",
    def=function(object) {standardGeneric("getNames_stacks")}
)

setMethod(
    f="getNames_stacks", 
    signature="Composante",
    definition=function(object) {
        return(object@names_stacks)
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
        if(getNameComp(object) == "niche_k" || getNameComp(object) == "niche_r" || getNameComp(object) == "generation") {
            flag = -1
            while(flag == -1) {
                cat("[Type 0 to quit] Which layer do you want to add?\n")
                cat(paste(1:length(object@remain_layers), ":", object@remain_layers))
                scanner = as.integer(readline())
                if(is.na(scanner) || scanner>length(object@remain_layers) || scanner<0){
                    print("ERROR: Your entry is incorrect, please try again")
                } else if(scanner == 0){
                    stop("You have stopped the program")
                }else{
                    object@nbModel = getNbModel(object)+1
                    newMod = model(getNameComp(object), getNbModel(object), getMethodComp(object), object@remain_layers[scanner])
                    object@listModel = c(object@listModel, newMod)
                    flag = 1
                }
            }
            
            
        } else {
            for(i in 1:nbToAdd) {
                object@nbModel = getNbModel(object)+1
                newMod = model(getNameComp(object), getNbModel(object), getMethodComp(object), "STANDARD")
                object@listModel = c(object@listModel, newMod)
            }
        }
        rm(newMod)
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
        # compteur = 0
        # numModelToDel = sort(numModelToDel)
        # Deletion of the models
        # for(i in numModelToDel) {
            # object@remain_layers = c(object@remain_layers,getNameModel(object@listModel[[i-compteur]])[3])
            # object@listModel = object@listModel[-(i-compteur)]
            # compteur = compteur+1
        # }
        object@remain_layers = c(object@remain_layers,getNameModel(object@listModel[[numModelToDel]])[3])
        object@listModel = object@listModel[-numModelToDel]
        object@nbModel = getNbModel(object) - length(numModelToDel)
        # Update of the models number
        for(i in 1: getNbModel(object) ) {
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
            if(getNameComp(object) == "niche_k" || getNameComp(object) == "niche_r" || getNameComp(object) == "generation") {
                cat("[Type 0 to quit] What do you want to do?\n1: Add a model\n2: Delete a model\n3: Change a model\n4: Change the independant model")
            } else {
                cat("[Type 0 to quit] What do you want to do?\n1: Add a model\n2: Delete a model\n3: Change a model")
            }
            scanner = as.integer(readline())
            if(is.na(scanner) || scanner>=5 || scanner<0){
                print("ERROR: Your entry is incorrect, please try again")
            } else if(scanner == 0){
                stop("You have stopped the program")
            }else{
                if(scanner == 4 && getNameComp(object) != "niche_k" && getNameComp(object) != "niche_r" && getNameComp(object) != "generation") {
                    print("ERROR: Your entry is incorrect, please try again")
                } else {
                    flag = 1
                }
            }
        }
        # Add a model
        if(scanner == 1){
            flag = 0
            while(flag == 0){
                if(getNameComp(object) == "niche_k" || getNameComp(object) == "niche_r" || getNameComp(object) == "generation") {
                    if(getNb_stacks(object) == getNbModel(object)) {
                        cat("ERROR: You can not have more models than stacks.")
                        flag = 1
                    } else {
                        cat("How many models do you want to add ?")
                        nbToAdd = as.integer(readline())
                        if(is.na(nbToAdd) || nbToAdd<1){
                            cat("ERROR: Your entry is incorrect, please try again.\n")
                        } else if(getNb_stacks(object) < getNbModel(object)+nbToAdd) {
                            cat("ERROR: You can not have more models than stacks.\n")
                        } else {    
                            flag = 1
                            object = addModel(object, nbToAdd)
                        }
                    }
                } else {
                    cat("How many models do you want to add ?")
                    nbToAdd = as.integer(readline())
                    if(is.na(nbToAdd) || nbToAdd<1){
                        print("ERROR: Your entry is incorrect, please try again.\n")
                    } else {    
                        flag = 1
                        object = addModel(object, nbToAdd)
                    }
                }
            }
        # Delete a model
        } else if(scanner == 2){
            flag = 0
            while(flag == 0){
                if(getNbModel(object) > 1) {
                    print(object)
                    cat("Which model do you want to delete ? (You can not delete the independant model)")
                    nbToDel = as.integer(readline())
                    if(is.na(nbToDel) || nbToDel> getNbModel(object) || nbToDel<1){
                        print("ERROR: Your entry is incorrect, please try again")
                    }else{
                        flag = 1
                        object = delModel(object, nbToDel)
                    }
                } else {
                    stop("There is only one model left in this composante. You can not delete it.")
                }
            }
        # Change a model
        } else if(scanner == 3) {
            flag = 0
            while(flag == 0){
                print(object)
                cat("Which model do you want to change (except for the independant model) ? Please enter the model's number.")
                change = as.integer(readline())
                if(is.na(change) || change> getNbModel(object) || change<1){
                    print("ERROR: Your entry is incorrect, please try again")
                }else{
                    flag = 1
                }
            }
            flag = 0
            while(flag == 0){
                cat("What do you want to change?\n1: Model function \n2: Model parameters")
                choice = as.integer(readline())
                if(is.na(choice) || choice> 2 || choice<1){
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
        # Change the independant model
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
        cat(getType_combinaison(object),"\n")
        if(getNameComp(object) == "niche_r" || getNameComp(object) == "niche_k" || getNameComp(object) == "generation") {
            cat(".................... Independant model ....................\n")
            cat("           ")
            print(object@independance)
        }
        for(i in 1:length(object@listModel)) {
            if(getNameComp(object) == "niche_r" || getNameComp(object) == "niche_k" || getNameComp(object) == "generation") {
                cat("............... Model:",i,",",object@listModel[[i]]@type_model,", name:",object@listModel[[i]]@name[3],"...............\n")
                print(object@listModel[[i]])
            } else {
                cat("............... Model:",i,",",object@listModel[[i]]@type_model,"...............\n")
                print(object@listModel[[i]])
            }
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
        if(getNameComp(object) == "niche_k" || getNameComp(object) == "niche_r" || getNameComp(object) == "generation") {
            if(all == 0) { 
                # ask which model the user wants to assess
                cat("[Type 0 to quit] Which model do you want to assess? Type the model number.\n")
                print(object)
                cat("Type", getNbModel(object)+1, "to assess the independant model.\n")
                cat("Type", getNbModel(object)+2, "to assess all models.\n")
                

                flag = -1
                while(flag == -1) {
                    choice = as.integer(readline())
                    if(choice!=0 && choice<=getNbModel(object)+2 && choice>0 && !is.na(choice)) {
                        if(choice == getNbModel(object)+1) {
                            object@independance = setResult_prior(object@independance)
                        } else if(choice == getNbModel(object)+2) {
                            object@independance = setResult_prior(object@independance)
                            for(i in 1:getNbModel(object)) {
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
                for(i in 1:getNbModel(object)) {
                    object@listModel[[i]] = setResultPriorMod(object@listModel[[i]],1)
                }
            }

        # if there is no independant model            
        } else {
            if(all == 0) { 
                # ask which model the user wants to assess
                cat("[Type 0 to quit] Which model do you want to assess? Type the model number.\n")
                print(object)
                cat("Type", getNbModel(object)+1, "to assess all models.\n")
                

                flag = -1
                while(flag == -1) {
                    choice = as.integer(readline())
                    if(choice!=0 && choice<=getNbModel(object)+1 && choice>0 && !is.na(choice)) {
                        if(choice == getNbModel(object)+1) {
                            for(i in 1:getNbModel(object)) {
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
                for(i in 1:getNbModel(object)) {
                    object@listModel[[i]] = setResultPriorMod(object@listModel[[i]],1)
                }
            }
        }
        return(object)
    }
)