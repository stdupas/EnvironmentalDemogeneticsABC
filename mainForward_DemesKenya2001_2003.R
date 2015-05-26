### Setting working directories
#rm(list=ls())
#wd="/media/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # portable
#wd="/media/dupas/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # fixe
#wd="/home/legs/GraphPOP/" # portable steph
#wd="/home/dupas/GraphPOP/" # fixe
#wd="/home/arno/These/GraphPOP" # portable arno
#wd="/home/arnaudb/Documents/GraphPOP" # labo arno
#setwd(wd)

### Sourcing functions files
source("Class_paramList.R")
source("Class_backward.R")
source("Class_forward.R")
source("PriorFunctions.R")
source("AskModelsFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("PriorFunctions.R")
source("forwardModels.R")
source("readNetCDF.R")
source("function_computeMeanEnvData.R")
source("ParallelGibbs.R")
source("script_distance_km_test.R")

### Sourcing Libraries
library(raster)
library(ape)
library(stringr)
library(lattice)
library(parallel)
library(lubridate) # of little use, should be removed (function "days()")
library(RNetCDF) # to read NetCDF data (climate spatial time series)
library(argosfilter)
library(Matrix)
library(foreach)
library(doParallel)

########### Parameters initialisation  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#Lecture et aggregation des fichiers de données
#EnvDataRasterStack = nc2EnvDataAndRasterStack(ncDirectory="../ForwardSimulData/",aggregationParam=40)
#Lecture des données déjà aggrégées. Gain de temps
EnvDataRasterStack = readRDS("../ForwardSimulData/ObjectEnvdataRasterStackAggr40")
rasterStack <- EnvDataRasterStack[[2]]
EnvData <- EnvDataRasterStack[[1]]
rm(EnvDataRasterStack)

###### release data 
release=as.data.frame(xyFromCell(rasterStack,1:ncell(rasterStack)))
birthDates <- as.Date(as.Date("2000/01/01"):as.Date("2000/06/01"),origin="1970/01/01")
release <- release[rep(1:nrow(release),length(birthDates)),]
release$birthDate=rep(birthDates,each=ncell(rasterStack))
release$size=1
release$demeNb <-  cellFromXY(object = rasterStack, xy = release[, c("x", "y")])
#release <- aggregation(release,BY=c("birthDate","demeNb"),methodes=c("Mean","Mean","Name","Sum","Name"))

###### recovery data 
recovery <- read.csv("../ForwardSimulData/Stemborer_Kenya2001_2005.csv", header = TRUE, sep = ";")
recovery <- recovery[,c("Long_dec","Lat_dec","Diss_Date","B._fusca", "no._plants")]
recovery = recovery[-which(is.na(recovery$no._plants)),]
recovery$B._fusca = as.numeric(recovery$B._fusca)/recovery$no._plants
recovery = recovery[,c("Long_dec","Lat_dec","Diss_Date","B._fusca")]
colnames(recovery) <- c("x","y","birthDate","size")
recovery$birthDate <- as.Date(as.character(recovery$birthDate),format="%d/%m/%y")
recovery$demeNb <- cellFromXY(object = rasterStack, xy = recovery[, c("x", "y")])
recovery <- aggregation(recovery,BY=c("birthDate","demeNb"),methodes=c("Mean","Mean","Name","Sum","Name"))
recovery <- recovery[which(as.Date(recovery$birthDate)<=as.Date("2003/12/31")),]

minDates = min(release$birthDate)
maxDates = min(max(recovery$birthDate),max(as.Date(colnames(EnvData))))
Dates <- as.Date(as.Date(minDates):as.Date(maxDates),origin="1970/01/01")

#Calcul de la distance en km entre deux demes à partir de leurs coordonées
distMat <- distanceMatrixFromRaster2(object = rasterStack)

demeSizes <- array(0,dim=c(nrow(EnvData),length(Dates)),dimnames = list(1:nrow(EnvData),as.character(Dates)))
demeSizes[,as.character(birthDates)] <- 10
#On récupère les dates d'interêt seulement.
EnvData <- EnvData[,colnames(demeSizes),]
#On fait la moyenne des precipitations sur les 15 derniers jours
EnvData2 = computeMeanEnvData(EnvData, "pr", 15)

#!!!!!!!!!!!!!!!!!!!!!!   POUR LES DONNEES REELLES, S'ARRETER ICI  !!!!!!!!!!!!!!!!!!!!!!!!!

#Formatage pour données de test:
#******************************

#Aggregation des jours par 10
EnvData2 = aggregateDays(EnvData2,10,Dates)
#On ne garde qu'une date sur 10
Dates = Dates[seq(10,length(Dates),by=10)]
#On ne garde les birthdates qui sont présentes dans Dates
tmp = birthDates %in% Dates
birthDates = birthDates[which(tmp == TRUE)]

#Création des données de test en lancant la fonction expectedInd sur des paramètres imposés.
recovery2 = buildDataSet()




#Fonctions de visualisation:
#***************************

plotParallelGibbs <- function(a=a,obj=1) {
    # obj=1 pour plot les parametres, obj=2 pour plot les posteriors
    if(obj==1) {
        par(mfrow=c(2,2))
        plot(a[[obj]][,1],t="l",ylab="K.pr.Xmin", xlab="Iteration",col="orange",ylim=c(0,2))
        points(a[[obj+2]][,1],t="l",col="navyblue")
        
        plot(a[[obj]][,5],t="l",ylab="R.pr.Xmin", xlab="Iteration",col="orange",ylim=c(0,2))
        points(a[[obj+2]][,5],t="l",col="navyblue")        
        
        plot(a[[obj]][,2],t="l",ylab="K.pr.Xmax", xlab="Iteration",col="orange",ylim=c(5,15))
        points(a[[obj+2]][,2],t="l",col="navyblue")
        
        plot(a[[obj]][,6],t="l",ylab="R.pr.Xmax", xlab="Iteration",col="orange",ylim=c(5,15))
        points(a[[obj+2]][,6],t="l",col="navyblue")
        
        plot(a[[obj]][,3],t="l",ylab="K.pr.Xopt", xlab="Iteration",col="orange",ylim=c(2,8))
        points(a[[obj+2]][,3],t="l",col="navyblue")
        
        plot(a[[obj]][,7],t="l",ylab="R.pr.Xopt", xlab="Iteration",col="orange",ylim=c(2,8))
        points(a[[obj+2]][,7],t="l",col="navyblue")
        
        plot(a[[obj]][,4],t="l",ylab="K.pr.Yopt", xlab="Iteration",col="orange",ylim=c(15,25))
        points(a[[obj+2]][,4],t="l",col="navyblue")
        
        plot(a[[obj]][,8],t="l",ylab="R.pr.Yopt", xlab="Iteration",col="orange",ylim=c(5,15))
        points(a[[obj+2]][,8],t="l",col="navyblue")      
        
        plot(a[[obj]][,9],t="l",ylab="R.tas.Xmin", xlab="Iteration",col="orange",ylim=c(270,295))
        points(a[[obj+2]][,9],t="l",col="navyblue")
        
        plot(a[[obj]][,10],t="l",ylab="R.tas.Xmax", xlab="Iteration",col="orange",ylim=c(295,315))
        points(a[[obj+2]][,10],t="l",col="navyblue")
        
        plot(a[[obj]][,11],t="l",ylab="R.tas.Xopt", xlab="Iteration",col="orange",ylim=c(285,305))
        points(a[[obj+2]][,11],t="l",col="navyblue")
        
        plot(a[[obj]][,12],t="l",ylab="R.tas.Yopt", xlab="Iteration",col="orange",ylim=c(0,2))
        points(a[[obj+2]][,12],t="l",col="navyblue") 
        
        
    } else {
        par(mfrow=c(1,1))
        plot(a[[obj]][,1],t="l", ylab="posteriors", xlab="Iteration", col="orange", ylim=c(-1300,-700))
        points(a[[obj+2]][,1], t="l", col="navyblue")
    }
}

plotAllValues <- function(a=a,burning=0) {
        obj=1

        intervalle = c(1:burning)

        par(mfrow=c(2,2))
        plot(a[[obj+4]][,1], t="l",ylab="K.pr.Xmin", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,1],prob=TRUE,xlab="K.pr.Xmin",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,1]),col="red")

        plot(a[[obj+4]][,5],t="l",ylab="R.pr.Xmin", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,5],prob=TRUE,xlab="R.pr.Xmin",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,5]),col="red")

        plot(a[[obj+4]][,2],t="l",ylab="K.pr.Xmax", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,2],prob=TRUE,xlab="K.pr.Xmax",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,2]),col="red")
        
        plot(a[[obj+4]][,6],t="l",ylab="R.pr.Xmax", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,6],prob=TRUE,xlab="R.pr.Xmax",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,6]),col="red")
        
        plot(a[[obj+4]][,3],t="l",ylab="K.pr.Xopt", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,3],prob=TRUE,xlab="K.pr.Xopt",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,3]),col="red")
        
        plot(a[[obj+4]][,7],t="l",ylab="R.pr.Xopt", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,7],prob=TRUE,xlab="R.pr.Xopt",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,7]),col="red")
        
        plot(a[[obj+4]][,4],t="l",ylab="K.pr.Yopt", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,4],prob=TRUE,xlab="K.pr.Yopt",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,4]),col="red")
        
        plot(a[[obj+4]][,8],t="l",ylab="R.pr.Yopt", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,8],prob=TRUE,xlab="R.pr.Yopt",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,8]),col="red")
        
        plot(a[[obj+4]][,9],t="l",ylab="R.tas.Xmin", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,9],prob=TRUE,xlab="R.tas.Xmin",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,9]),col="red")
        
        plot(a[[obj+4]][,10],t="l",ylab="R.tas.Xmax", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,10],prob=TRUE,xlab="R.tas.Xmax",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,10]),col="red")
        
        plot(a[[obj+4]][,11],t="l",ylab="R.tas.Xopt", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,11],prob=TRUE,xlab="R.tas.Xopt",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,11]),col="red")
        
        plot(a[[obj+4]][,12],t="l",ylab="R.tas.Yopt", xlab="Iteration",col="orange")
        hist(a[[obj+4]][-intervalle,12],prob=TRUE,xlab="R.tas.Yopt",breaks=50,main="")
        lines(density(a[[obj+4]][-intervalle,12]),col="red")

}
