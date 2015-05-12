aggregation <- function(df,BY,methodes)
{
    # Function to aggregate lines of a data.frame
    # arguments
    # df: data.frame to aggregate
    # BY: columns to use for aggregation of lines, lines having the same value in all the columns will be aggregated
    # value
    # a data.frame aggregated. Numeric and integer columns are given mean value of aggregated lines. Character and other class
    # of vectors are given collapsed character string value, or first value if all the same
    
    df$MergeBY <- NA
    for (ligne in 1:nrow(df)){df$MergeBY[ligne]=paste(df[ligne,BY],collapse="_")}
    df <- df[order(df[,"MergeBY"]),];rownames(df)=1:nrow(df)
    df_return=df[1,];df_return=df_return[-1,]
    ligne=1
    methodes = append(methodes,"Name")
    while (ligne <= dim(df)[1])
    {
        print(ligne)
        lignes = which(df[,"MergeBY"]==df[ligne,"MergeBY"])
        if (length(lignes)==0) {ligne=ligne+1} else{
            sub_df <- df[lignes,]
            if (length(lignes)>1) {df_return <- rbind(df_return,sub_df[1,])
                                   for (colonne in 1:dim(sub_df)[2]) {
                                       df_return[nrow(df_return),colonne] = switch(methodes[colonne],
                                                                                   Mean = mean(na.omit(sub_df[,colonne])),
                                                                                   Sum =  sum(na.omit(sub_df[,colonne])),
                                                                                   Name = na.omit(sub_df[,colonne])[1],
                                                                                   Paste = paste(sub_df[,colonne],collapse="_")
                                       )
                                   } 
            }
            else {df_return <- rbind(df_return,df[lignes,])}
            ligne=ligne+length(lignes)
        }
    }
    return(df_return[,-ncol(df_return)])
}

forwardMigrationRateMatrixFromKernel <- function(dispersionKernel){
    # Normalizes a matrix of dispersion kernel between cells so that individuals bounce on the frontier.
    #
    # Args:
    #   dispersion: a matrix representing the values of a specified kernel (function of distances between cells)
    #
    # Returns:
    #   A migration rate matrix (note that rowSums and colSums are not 1: cause of bordure effect, individuals go "out of the world")
    # ex: 
    # kernelMatrix=c(1.000000e+00,0.804259048,0.004416269,8.156667e-05,
    #                 8.042590e-01,1.000000000,0.804259048,4.416269e-03,
    #                 4.416269e-03,0.804259048,1.000000000,8.042590e-01,
    #                 8.156667e-05,0.004416269,0.804259048,1.000000e+00),nrow=4,ncol=4)
    # fmatrix <- forwardMigrationRateMatrixFromKernel(kernelMatrix)
    # rowSums(fmatrix)
    return(dispersionKernel/matrix(rowSums(dispersionKernel),nrow=nrow(dispersionKernel),ncol=ncol(dispersionKernel),byrow=FALSE))
}

reproduction <- function(tipDemes, r)
{
    # Reproduces individuals according to growth rate of their deme
    # Arguments:
    # - tipeDemes : a data frame describing the individuals to reproduce
    # with columns "individualNb" (individual number) and "demeNb" (deme number)
    # - r : a vector of growth rate of each deme
    # Value:
    # A data.frame describing the offsprings, 
    # - column "individualNb" gives the number of the parent, 
    # - column "demeNb" gives the deme 
    # rownames have the format "parent_number.offspring_number"
    # 
    # Example
    # r<-2^(0:3)
    # tipDemes=data.frame(individualNb=1:10,demeNb=sample(1:4,10,replace=TRUE))
    # reproduction(tipDemes,r)
    rownames(tipDemes)<-tipDemes$individualNb
    new_individuals <- rep(tipDemes$individualNb,rpois(nrow(tipDemes),r[tipDemes$demeNb ]))
    newTipDemes <- tipDemes[new_individuals,]
    return(newTipDemes)
}

reproductionMigration <- function(r,migrationMatrix)
{
    # Calculates number of offspring in deme j of one individual in deme i
    #
    # Arguments:
    # - tipDemes : a data frame describing the individuals to reproduce
    # with columns "individualNb" (individual number) and "demeNb" (deme number)
    # - r : a vector of growth rate of each deme
    # - K : a vector of carrying capacity of each deme
    # Value:
    # A transiton matrix describing parent offspring reproduciton and movements
    # 
    # Example
    # r<-2^(0:3)
    # migrationMatrix <- matrix(c(.9,.1,0,0,.05,.9,0.05,0,0,.05,.9,0.05,0,0,.1,.9),nrow=4,ncol=4,byrow=TRUE)
    # reproductionMigration(r,K,migrationMatrix)
    #
    #
    
    r%*%migrationMatrix
}

individualsToRemoveFromCompet <- function(individuals,demeNb,K)
{
    # Samples individuals exceeding carring capacity for removal in a vector of deme attribution
    # 
    # CarryingCap=rep(2,4)
    # demeNumb=c(1,1,1,2,3,3,3,3,4,4,4)
    # indiv = 1:11
    # individualsToRemoveFromCompet(individuals=indiv,demeNb=demeNumb,K=CarryingCap)
    numberPerDeme <- hist(demeNb,plot=FALSE,breaks=(0:length(K))+.5)$counts
    excedent <- round((numberPerDeme - K > 0)*(numberPerDeme - K))
    remove = NULL
    for (i in 1:length(excedent))
    {
        remove=append(remove,sample(x=individuals[which(demeNb==i)],size=excedent[i]))
    }
    return(remove)
}

reprMigr <- function(tipDemes,migrationMatrix,r,generationTime,generationTimeRelativeSD)
{
    # Randomly reproduces and migrates individuals according to growth rate of their deme
    # migration rates to descendent demes and carrying capcacity of descendent dames
    # Arguments:
    # - tipeDemes : a data frame describing the individuals to reproduce
    # with columns "individualNb" (individual number), "demeNb" (deme number) and "birthDate"
    # - migrationMatrix : absolute transition matrix with expected offspring size in column for each deme in line
    # - r : reproduction rate in each deme
    # - generationTime : a vector of generation time for each deme
    # - generationTimeRelativeSD : a standard error relative to the mean for generation time
    # Value:
    # A data.frame describing the offsprings, 
    # - column "individualNb" gives the number of the parent, 
    # - column "demeNb" gives the deme 
    # - column "birthdate"
    # rownames have the format "parent_number.offspring_number"
    # 
    # Example
    # r<-2^(0:3)
    # tipDemes=data.frame(individualNb=1:10,birthDate=as.Date("2013/01/23"),demeNb=sample(1:4,10,replace=TRUE))
    # generationTime=sample(20:25,4,replace=TRUE)
    # generationTimeRelativeSD=.1
    # migrationMatrix <- matrix(c(.9,.1,0,0,.05,.9,0.05,0,0,.05,.9,0.05,0,0,.1,.9),nrow=4,ncol=4,byrow=TRUE) 
    # reprMigr(tipDemes=subset(tipDemes,tipDemes$birthDate=="2013-01-23"),
    #          migrationMatrix=migrationmatrix,
    #          r=r,
    #          generationTime,
    #          generationTimeRelativeSD)
    #
    
    #rownames(tipDemes)<-tipDemes$individualNb
    repMig = r*migrationMatrix
    # caculating absolute reproduction migration matrix with expected 
    # offspring size in column for each deme of the tip in line
    #trTip <- subset(repMig,subset=(1:nrow(repMig)==tipDemes$demeNb))
    trTip <- repMig[tipDemes$demeNb,]; if(class(trTip)=="numeric") {trTip <- as.matrix(t(trTip))}
    # generating number of offspring of each tip (line) in each deme (column)
    individualNb <- matrix( rpois(prod(dim(trTip)),trTip),nrow=nrow(trTip),ncol=ncol(trTip))
    # 
    offsTipDemes <- tipDemes[as.character(rep(rownames(tipDemes),rowSums(individualNb))),]
    #offsTipDemes$birthDate <- offsTipDemes$birthDate + rnorm(nrow(offsTipDemes),generationTime[offsTipDemes$demeNb],generationTimeRelativeSD*generationTime[offsTipDemes$demeNb])
    offsTipDemes$birthDate <- as.character(offsTipDemes$birthDate + rnorm(nrow(offsTipDemes),generationTime,generationTimeRelativeSD*generationTime))
    offsTipDemes$demeNb <- rep(rep(1:length(r),nrow(individualNb)),as.vector(t(individualNb)))
    return(offsTipDemes)
}

generate_parameterSeries <- function(EnvData,migrationMatrix,nicheKModelsList,nicheRModelsList,paramKList,paramRList)
{
    r <- K <- EnvData[,,1]
    for (Date in colnames(r)) #Date=colnames(r)[3]
    {
        values(rasterStack) <- EnvData[,Date,]
        # launch the siulation
        # Get the carrying capacity map :
        r[,Date] <- values(nicheFunctionForRasterStack(functionList = nicheKModelsList, 
                                                       rasterStack = rasterStack,
                                                       args = paramKList))
        
        # Get growth rate map :
        K[,Date] <- values(nicheFunctionForRasterStack(functionList = nicheRModelsList, 
                                                       rasterStack = rasterStack,
                                                       args = paramKList))
        
    }
    list(r=r,K=K)
}

generateIndividualSerie <- function(release,r_serie,K_serie,migrationMatrix, startingDate, stoppingDate,generationTime,generationTimeSD)
{
    # function to generate a time serie of individuals birth date and deme
    # arguments:
    # r_serie : matrix of growh rate per deme in line and day in column
    # K_serie : matrix of carrying capacity per deme in line and day in column
    # release : a data frame describing the individuals to reproduce
    # startingDate, stoppingDate: dates to start and stop the simulation
    # stoppingDate: date to stop the time serie simulation
    # with columns "individualNb" (individual number), "demeNb" (deme number) and "birthDate"
    # value:
    # a data frame with all individuals characteristics in the time serie
    # example:
    #
    # r_serie = matrix(rpois(731*4,c(3,4,5,7)),
    #                  nrow = 4, ncol = 731, byrow=TRUE,
    #                  dimnames = list(1:4,as.character(as.Date(as.Date("2001/01/01"):as.Date("2003/01/01"),origin="1970-01-01"))))
    # K_serie = matrix(10,nrow = 4, ncol = 731, byrow=TRUE,
    #                  dimnames = list(1:4,as.character(as.Date(as.Date("2001/01/01"):as.Date("2003/01/01"),origin="1970-01-01"))))
    # release = data.frame(individualNb=1:10,demeNb=sample(1:4,10,replace=TRUE),
    #                     birthDate=as.Date(sample(as.Date(as.Date("2001/01/01"):as.Date("2001/02/01"),
    #                                       origin="1970-01-01"),10,replace=TRUE)))
    # startingDate=as.Date("2001/01/01")
    # stoppingDate=as.Date("2003/01/01")
    # generationTime=21,generationTimeSD=3
    individuals <- release
    for (Date in as.Date(startingDate:stoppingDate,origin="1970-01-01"))
    {
        absoluteForwardTransitionMatrix <- absoluteForwardTransition(r_serie[,as.character(Date)],K_serie[,as.character(Date)],migrationMatrix)
        individuals <- rbind(individuals,reprMigrCompet(individuals[which(individuals$birthDate==Date),],absoluteForwardTransitionMatrix,generationTime,generationTimeSD))
    }
    individuals
}

demeReprodMigr <- function()
{
    
}

generateDemeSizeSerie <- function(release,
                                  rasterStack,
                                  dispersionFunction,
                                  dispersionParameters,
                                  nicheKFunctionList,
                                  nicheRFunctionList,
                                  nicheKParametersList,
                                  nicheRParametersList,
                                  generationTimeParameters,
                                  stoppingDate)
{
    kernelMatrix <- dispersionFunctionForRasterLayer(dispersionFunction=dispersionFunction,
                                                     rasterLayer=rasterStack[[1]],
                                                     args=dispersionParameters)
    migrationMatrix <- forwardMigrationRateMatrixFromKernel(kernelMatrix)
    rm(kernelMatrix)
    generationTime <- generationTimeParameters$mean
    generationTimeSD <- generationTimeParameters$mean*generationTimeParameters$SD
    release$demeNb <- cellFromXY(object = rasterStack, xy = release[, c("x", "y")])
    release$demeSize <- 1
    Dates = as.Date(min(release$birthDate):as.Date(stoppingDate),origin="1970-01-01")
    demeSizes <- array(0,dim=c(nrow(EnvData),length(Dates)),dimnames = list(1:nrow(EnvData),as.character(Dates)))
    release <- aggregation(release,BY=c("birthDate","demeNb"),methodes=c("Paste","Mean","Mean","Name","Name","Sum"))
    for (i in rownames(release))
    {demeSizes[release[i,"demeNb"],as.character(release[i,"birthDate"])] <- release[i,"demeSize"]}
    for (Date in as.character(as.Date(as.Date(min(Dates)):as.Date(stoppingDate),origin="1970-01-01"))) #Date=as.character(as.Date(as.Date(min(Dates)):as.Date(fwParamList$stoppingDate),origin="1970-01-01"))[1]
    {
        values(rasterStack) <- EnvData[,as.character(Date),]
        K <- values(nicheFunctionForRasterStack(functionList = nicheKFunctionList, 
                                                rasterStack = rasterStack,
                                                args = nicheKParametersList
        )
        )
        # competition: cuts deme sizes to K
        demeSizes[,Date] <- (demeSizes[,as.character(Date)]>K)*K + round((demeSizes[,as.character(Date)]<=K)*demeSizes[,as.character(Date)])
        # reproduction growth rate
        r <- values(nicheFunctionForRasterStack(functionList = nicheRFunctionList, 
                                                rasterStack = rasterStack,
                                                args = nicheRParametersList
        )
        )
        descendants <- rpois(nrow(demeSizes),demeSizes[,as.character(Date)]*r)
        #migration
        descendants <- rpois(nrow(demeSizes),descendants%*%migrationMatrix)
        descendantDates <- as.character(as.Date(Date) + rnorm(sum(descendants),generationTime,generationTimeSD))
        demes <- rep(1:nrow(demeSizes),descendants)
        if (length(descendantDates>0)) {
            for (i in 1:length(descendantDates))
            {
                demeSizes[demes[i],as.character(descendantDates[i])] <- demeSizes[demes[i],descendantDates[i]] + 1
            }
        }
    }
    demeSizes  
}

demeSizeArray <- function(release, minDate, maxDate, BY="day")
{
    #
    # construction of expected deme sizes array
    #
    release$demeSize <- 1
    Dates = as.Date(as.Date(minDate):as.Date(maxDate),origin="1970-01-01")
    demeSizes <- array(0,dim=c(nrow(EnvData),length(Dates)),dimnames = list(1:nrow(EnvData),as.character(Dates)))
    release <- aggregation(release,BY=c("birthDate","demeNb"),methodes=c("Mean","Mean","Sum","Name","Name","Sum"))
    for (i in rownames(release))
    {demeSizes[release[i,"demeNb"],as.character(release[i,"birthDate"])] <- release[i,"demeSize"]}
    demeSizes  
}

table2arrayFromArrayModel <- function(Table=release,arrayModel=EnvData,Dates)
{ 
    #
    # construction of expected recovery array in the same grid as EnvData with subset of dates
    # Args:
    # Table : a table with columns "demeNb" corresponding to EnvData line
    #                              "size" corresponding to number of individual observed
    #                              "birthDate"corrsponding to date of observation
    #
    demeSizes <- array(0,dim=c(nrow(EnvData),length(Dates)),dimnames = list(1:nrow(EnvData),as.character(Dates)))
    for (i in rownames(Table))
    {demeSizes[recovery[i,"demeNb"],as.character(Table[i,"birthDate"])] <- Table[i,"size"]}
    demeSizes
}



likelihood <- function(dispersionParameters = list(alpha=50,beta=2),
                       nicheKParametersList=list(pr=list(X0=0,Xopt=267.1267,Yopt=5),
                                                 tasmax=list(X0=273.1523,Xopt=321.9696,Yopt=1),
                                                 tasmin=list(X0=268.2348,Xopt=304.5951,Yopt=1)),
                       nicheRParametersList=list(pr=list(X0=0,Xopt=267.1267,Yopt=20),
                                                 tasmax=list(X0=273.1523,Xopt=321.9696,Yopt=1),
                                                 tasmin=list(X0=268.2348,Xopt=304.5951,Yopt=1)),
                       generationTimeParameters=list(mean=25,SD=3)
)
{
    
    #
    # building migration matrix
    #
    
    dispersionKernel <- fatTail1(distMat, dispersionParameters$alpha, dispersionParameters$beta)
    migrationMatrix <- migrationRateMatrix(dispersionKernel)  
    
    rm(dispersionKernel)
    
    #
    # Probability density of generation time inthe interval [mean-3SD,mean+3SD]
    #
    
    generationTimeInterval <- (generationTimeParameters$mean-3*round(generationTimeParameters$SD)):(generationTimeParameters$mean+3*round(generationTimeParameters$SD))
    generationTimeDensity <- dnorm(generationTimeInterval,generationTimeParameters$mean,generationTimeParameters$SD)
    generationTimeReducedInterval <- 1:length(generationTimeInterval)
    
    # construction of likelihood with expected recovery
    for (i in 1:(ncol(demeSizes)-max(generationTimeInterval))) # Date = colnames(demeSizes)[1]
    {
        values(rasterStack) <- EnvData[,i,]
        K <- values(nicheFunctionForRasterStack(functionList = nicheKFunctionList, 
                                                rasterStack = rasterStack,
                                                args = nicheKParametersList,
                                                meth = "product")
        )
        R <- values(nicheFunctionForRasterStack(functionList = nicheRFunctionList, 
                                                rasterStack = rasterStack,
                                                args = nicheRParametersList,
                                                meth = "product")
        )
        R[is.na(R)]<-0
        K[is.na(K)]<-0
        # reproduction: multiplies by r
        reproducedAtDate <- demeSizes[,i]*R
        # competition: cuts deme sizes to K
        competedAtDate <- (reproducedAtDate>K)*K + (reproducedAtDate<=K)*reproducedAtDate
        # migration: moves as adult after development (after generation time interval)
        competedAtDate <- competedAtDate%*%migrationMatrix
        demeSizes[,i+generationTimeInterval] <- demeSizes[,i+generationTimeInterval]+t(competedAtDate)%*%t(generationTimeDensity)
    }
    logLikelihood <- sum(dpois(recovery[,"size"],demeSizes[recovery[,"demeNb"],as.character(recovery[,"birthDate"])],log=TRUE))
}


likelihoodShort <- function(dispersionRate = .025,dispersionDistance=100,
                            K.pr.X0=0,K.pr.Xopt=267.1267,K.pr.Yopt=5,K.generationTime=25,K.generationTimeSD=3,
                            R.pr.X0=0,R.pr.Xopt=267.1267,R.pr.Yopt=5,R.generationTime=25,R.generationTimeSD=3)
{
    
    #
    # building migration matrix
    #
    
    migrationMatrix <- (!(distMat == 0)&(distMat < dispersionDistance))*dispersionRate + (distMat==0)*(1-dispersionRate*4)
    migrationMatrix <- migrationMatrix/colSums(migrationMatrix)
    
    #
    # Probability density of generation time inthe interval [mean-3SD,mean+3SD]
    #
    
    generationTimeInterval <- (generationTime-3*round(generationTimeSD)):(generationTime+3*round(generationTimeSD))
    generationTimeDensity <- dnorm(generationTimeInterval,generationTime,generationTimeSD)
    generationTimeReducedInterval <- 1:length(generationTimeInterval)
    
    # construction of likelihood with expected recovery
    system.time(
        for (i in 1:(ncol(demeSizes)-max(generationTimeInterval))) # Date = colnames(demeSizes)[1]
        {
            K <- linearTreeParameters(EnvData[,i,"pr"],K.pr.X0,K.pr.Xopt,K.pr.Yopt)
            R <- linearTreeParameters(EnvData[,i,"pr"],R.pr.X0,R.pr.Xopt,R.pr.Yopt)
            R[is.na(R)]<-0
            K[is.na(K)]<-0
            # reproduction: multiplies by r
            reproducedAtDate <- demeSizes[,i]*R
            # competition: cuts deme sizes to K
            competedAtDate <- (reproducedAtDate>K)*K + (reproducedAtDate<=K)*reproducedAtDate
            # migration: moves as adult after development (after generation time interval)
            competedAtDate <- competedAtDate%*%migrationMatrix
            demeSizes[,i+generationTimeInterval] <- demeSizes[,i+generationTimeInterval]+t(competedAtDate)%*%t(generationTimeDensity)
        }
    )
    logLikelihood <- sum(dpois(recovery[,"size"],demeSizes[recovery[,"demeNb"],as.character(recovery[,"birthDate"])],log=TRUE))
    logLikelihood
}





#############################################################
#############################################################
##                                                         ##
##                     FONCTIONS BAYESIAN                  ##
##                                                         ##
#############################################################
#############################################################


logPostDens <- function(start){
    # Fonction calculant le posterior pour un lot de parametres donne
    # Variables: 
    #           start: vecteur contenant les valeurs des parametres
    # Retourne:
    #           somme du log de la likelihood et du log du prior

    # Affectation des vairables puis alcul du loglikelihood
    K.pr.Xmin = start[1]
    K.pr.Xmax = start[2]
    K.pr.Xopt = start[3]
    K.pr.Yopt = start[4]
    R.pr.Xmin = start[5]
    R.pr.Xmax = start[6]
    R.pr.Xopt = start[7]
    R.pr.Yopt = start[8]
    R.tas.Xmin = start[9]
    R.tas.Xmax = start[10]
    R.tas.Xopt = start[11]
    R.tas.Yopt = start[12]
    loglike = likelihoodShortTest(K.pr.Xmin, K.pr.Xmax, K.pr.Xopt, K.pr.Yopt,
                                  R.pr.Xmin, R.pr.Xmax, R.pr.Xopt, R.pr.Yopt,
                                  R.tas.Xmin, R.tas.Xmax, R.tas.Xopt, R.tas.Yopt)
    
    # Calcul des priors
    pKXmin = dunif(K.pr.Xmin, min=0, max=2) + 10^-320
    pKXmax = dunif(K.pr.Xmax, min=5, max=15) + 10^-320
    pKXopt = dunif(K.pr.Xopt, min=2, max=8) + 10^-320
    pKYopt = dunif(K.pr.Yopt, min=15, max=25) + 10^-320
    pRXmin = dunif(R.pr.Xmin, min=0, max=2) + 10^-320
    pRXmax  = dunif(R.pr.Xmax, min=5, max=15) + 10^-320
    pRXopt = dunif(R.pr.Xopt, min=2, max=8) + 10^-320
    pRYopt = dunif(R.pr.Yopt, min=5, max=15) + 10^-320
    tRXmin = dunif(R.tas.Xmin, min=270, max=295) + 10^-320
    tRXmax = dunif(R.tas.Xmax, min=296, max=315) + 10^-320
    tRXopt = dunif(R.tas.Xopt, min=285, max=305) + 10^-320
    tRYopt = dunif(R.tas.Yopt, min=0, max=2) + 10^-320
    logprior = sum(sapply(c(pKXmin,pKXmax,pKXopt,pKYopt,
                            pRXmin,pRXmax,pRXopt,pRYopt,
                            tRXmin,tRXmax,tRXopt,tRYopt),
                          FUN = log))
    
    return(loglike + logprior)
    
}


likelihoodShortTest <- function(
    K.pr.Xmin=0.5, K.pr.Xmax=10, K.pr.Xopt=4, K.pr.Yopt=20,
    R.pr.Xmin=0.5, R.pr.Xmax=10, R.pr.Xopt=4, R.pr.Yopt=10,
    R.tas.Xmin=285, R.tas.Xmax=305, R.tas.Xopt=295, R.tas.Yopt=1){
    # Fonction qui calcule la log likelihood pour pour un lot de parametres donne
    # Variables: 
    #           K.pr: parametres pour la fonction K=f(precipitation)
    #           R.pr: parametres pour la fonction R=f(precipitation)
    #           R.tas: parametres pour la fonction R=f(temperature)

    # Calcul du nombre d'individus attendus pour ce lot de parametres
    larveSizes = expectedInd(K.pr.Xmin, K.pr.Xmax, K.pr.Xopt, K.pr.Yopt,
                             R.pr.Xmin, R.pr.Xmax, R.pr.Xopt, R.pr.Yopt,
                             R.tas.Xmin, R.tas.Xmax, R.tas.Xopt, R.tas.Yopt)
    
    # On ne garde que les donnees avec un couple (deme,date) commun a notre recovery
    result = larveSizes[cbind(recovery2[,"demeNb"], as.character(recovery2[,"birthDate"]))]
    
    # Si recovery est > 0 et si result est egal à 0, dpois retourne -Inf
    # Il faut donc convertir les 0 de result en 0.0001 (ou autre different de 0)
    result[which(result == 0)] = 0.0001
    result[is.nan(result)] = 0.0001
    
    # Si recovery n'est pas un integer, dpois retourne -Inf
    # Il faut donc arrondir les valeurs de recovery
    logLikelihood <- sum(dpois(round(recovery2[,"size"]) , result , log=TRUE))
    return(logLikelihood)
}


expectedInd <- function(
    K.pr.Xmin=0.5, K.pr.Xmax=10, K.pr.Xopt=4, K.pr.Yopt=20,
    R.pr.Xmin=0.5, R.pr.Xmax=10, R.pr.Xopt=4, R.pr.Yopt=10,
    R.tas.Xmin=285, R.tas.Xmax=305, R.tas.Xopt=295, R.tas.Yopt=1){
    # Fonction qui calcule le nombre d'individus attendus
    # Variables: 
    #           K.pr: parametres pour la fonction K=f(precipitation)
    #           R.pr: parametres pour la fonction R=f(precipitation)
    #           R.tas: parametres pour la fonction R=f(temperature)

    
    # Affectation des valeurs pour la migration, 
    # le temps de developpement et le temps de generation
    dispersionRate = .025;dispersionDistance=300;      
    
    generationTime = ceiling(25/10);
    generationTimeSD=ceiling(3/10);    
    dvlpTime=1+ceiling(5/10);
    dvlpTimeSD=1;
    
    # generationTime = 25;
    # generationTimeSD = 3;    
    # dvlpTime= 5;
    # dvlpTimeSD=1;

    # Matrice contenant les individus à l'extérieur des mais.
    parentSizes <- array(0,dim=c(nrow(EnvData2),length(Dates)),dimnames = list(1:nrow(EnvData2),as.character(Dates)))
    parentSizes[,as.character(birthDates)] <- 1
    
    # Matrice des individus à l'intérieur des mais.
    larveSizes <- array(0,dim=c(nrow(EnvData2),length(Dates)),dimnames = list(1:nrow(EnvData2),as.character(Dates)))
    
    # Matrice de migration
    migrationMatrix = Matrix(0, nrow = dim(distMat)[1], ncol = dim(distMat)[2], sparse = TRUE)
    ind1 = which((distMat != 0) & (distMat < dispersionDistance))
    ind2 = which(distMat == 0)
    
    migrationMatrix[ind1] = dispersionRate
    migrationMatrix[ind2] = 1-dispersionRate*4
    
    # Densite de probabilite du temps de generation dans un intervalle
    generationTimeInterval <- (generationTime-generationTimeSD):(generationTime+generationTimeSD)
    generationTimeDensity <- dnorm(generationTimeInterval,generationTime,generationTimeSD)
    
    dvlpTimeInterval <- (dvlpTime-dvlpTimeSD):(dvlpTime+dvlpTimeSD)
    dvlpTimeDensity <- dnorm(dvlpTimeInterval,dvlpTime,dvlpTimeSD)
    
    # Calcul du nombre d'individus attendus en fonction des parametres
    for (i in 2:(ncol(larveSizes)-max(generationTimeInterval))) 
    {
        # Calcul des parametres K et R en fonction des parametres et de la precipitation + temperature
        R.tasmin <- conquadraticSkewed1(EnvData2[,i,"tasmin"], R.tas.Xmin, R.tas.Xmax, R.tas.Xopt, R.tas.Yopt)
        R.tasmax <- conquadraticSkewed1(EnvData2[,i,"tasmax"], R.tas.Xmin, R.tas.Xmax, R.tas.Xopt, R.tas.Yopt)
        R.pr <- conquadraticSkewed1(EnvData2[,i,"pr"], R.pr.Xmin, R.pr.Xmax, R.pr.Xopt, R.pr.Yopt)
        R <- R.tasmax*R.tasmin*R.pr
        
        K <- conquadraticSkewed1(EnvData2[,i,"pr"], K.pr.Xmin, K.pr.Xmax, K.pr.Xopt, K.pr.Yopt)  
        
        R[is.na(R)]<-0
        K[is.na(K)]<-0
        
        R[is.nan(R)]<-0
        K[is.nan(K)]<-0
        
        # Migration des adultes
        migratedAtDate = parentSizes[,(i-1)]%*%migrationMatrix
        parentSizes[,i] = parentSizes[,i] + migratedAtDate[1,]
        parentSizes[which(parentSizes[,i]<0),i] = 0
        
        # Reproduction des adultes
        nbNaissancesOld = parentSizes[,i]*R
        nbNaissances = nbNaissancesOld
        
        tmp = larveSizes[,i-1] + nbNaissancesOld + larveSizes[,i]
        ind = which(tmp >= K)
        nbNaissances[ ind ] = (K[ind] - larveSizes[ind,i-1] - larveSizes[ind,i])*((K[ind] - larveSizes[ind,i-1] - larveSizes[ind,i])>0)
        
        larveSizes[,i] = larveSizes[,i-1] + nbNaissances + larveSizes[,i]
        larveSizes[which(larveSizes[,i]<0),i] = 0
        
        # Programmation de leur eclosion en papillon
        larveSizes[,i+generationTimeInterval] =  larveSizes[,i+generationTimeInterval] - outer(nbNaissances,generationTimeDensity,"*")
        parentSizes[,i+generationTimeInterval] = parentSizes[,i+generationTimeInterval] + outer(nbNaissances,generationTimeDensity,"*")
        
        # Programmation de leur mort
        # Attention, suprression des adultes dans leur deme de naissance, ne prends pas en compte la migration, c'est pas bien...
        tempsVie = 2
        if(i+max(generationTimeInterval)+tempsVie <= dim(parentSizes)[2]){
            parentSizes[,i+tempsVie+generationTimeInterval] = parentSizes[,i+tempsVie+generationTimeInterval] - outer(nbNaissancesOld,generationTimeDensity,"*")
        } 
    }
    return(larveSizes)
}

buildDataSet <- function() {
    # Fonction qui simule un recovery (pour des tests de convergence)

    # Calcul des valeurs possibles avec la meme fonction que pour la likelihood
    possibleData = expectedInd(K.pr.Xmin=0.5, K.pr.Xmax=10, K.pr.Xopt=4, K.pr.Yopt=20,
                               R.pr.Xmin=0.5, R.pr.Xmax=10, R.pr.Xopt=4, R.pr.Yopt=10,
                               R.tas.Xmin=285, R.tas.Xmax=305, R.tas.Xopt=295, R.tas.Yopt=1)
    
    # Nombre de donnees a calculer
    nbData = 1200

    # On pioche certains demes et certaines dates
    choiceDate = sample(1:(length(Dates)-5), nbData, replace=TRUE, prob=NULL)
    choiceDeme = sample(1:dim(EnvData2)[1], nbData, replace=TRUE, prob=NULL)
    
    # On forme les couples (deme,date)
    pData = NULL
    for(i in 1:nbData) {
        pData = rbind(pData, possibleData[choiceDeme[i], choiceDate[i]]) 
    }
    
    # On calcule le nombre d'individus
    sizes = rpois(nbData, pData)
    dates = colnames(possibleData[,choiceDate])
    
    # On transforme les integers de dates en vraies dates
    dataSet = cbind.data.frame(as.Date(dates),as.numeric(sizes),as.integer(choiceDeme))
    colnames(dataSet) = c("birthDate", "size", "demeNb")
    
    return(dataSet)
}
