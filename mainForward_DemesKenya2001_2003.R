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
library("argosfilter")
library(Matrix)
library(foreach)
library(doParallel)

########### Parameters initialisation  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#EnvDataRasterStack = nc2EnvDataAndRasterStack(ncDirectory="/Users/Stagiaire/ForwardSimulData/",aggregationParam=40)
EnvDataRasterStack = readRDS("/Users/Stagiaire/ForwardSimulData/ObjectEnvdataRasterStack")
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
recovery <- read.csv("/Users/Stagiaire/ForwardSimulData/Stemborer_Kenya2001_2005.csv", header = TRUE, sep = ";")
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

dispersionModel= "fatTail1"

nicheKFunctionList=list(pr="linearTreeParameters",
                        tasmax="linearTreeParameters",
                        tasmin="linearTreeParameters")
nicheRFunctionList=list(pr="linearTreeParameters",
                        tasmax="linearTreeParameters",
                        tasmin="linearTreeParameters")

#recovery_sim <- data.frame(individualNb=1:3000,
#                       x=runif(3000,min=bbox(rasterStack)[1,1],
#                               max=bbox(rasterStack)[1,2]),
#                       y=runif(3000,min=bbox(rasterStack)[2,1],
#                               max=bbox(rasterStack)[2,2]),
#                       birthDate=as.Date("2001/01/01")+days(sample(100:730,3000,replace=TRUE)),
#                       size=sample(1:10,3000,replace=TRUE)
#                       )

distMat <- distanceMatrixFromRaster2(object = rasterStack)
demeSizes <- array(0,dim=c(nrow(EnvData),length(Dates)),dimnames = list(1:nrow(EnvData),as.character(Dates)))
demeSizes[,as.character(birthDates)] <- 10
EnvData <- EnvData[,colnames(demeSizes),]

EnvData2 = computeMeanEnvData(EnvData, "pr", 15)
EnvData2 = aggregateDays(EnvData2,10,Dates)
Dates = Dates[seq(10,length(Dates),by=10)]
tmp = birthDates %in% Dates
birthDates = birthDates[which(tmp == TRUE)]

recovery2 = buildDataSet()

LimiteLikelihood <- function(){
  ############# Boucle de visualisation fastidieuse ###############
  x = seq(33,37,0.1)
  y = seq(8,12,0.1)
  z = NULL
  for ( i in 1: length(x)){
    zbis = NULL
    cat("\n", i, ": ")
    for( j in 1: length(y)){
      cat("*")
      zbis = c(zbis, likelihoodShortTest(x[i], y[j]))
    }
    z = rbind(z, zbis)
  }
  zorder = rank(z)
  library(rgl)
  persp3d(x,y,z,col=rainbow(as.integer(max(zorder)))[zorder])
  filled.contour(x,y,z, color.palette = heat.colors)
  
}


plotGrosGibbs <- function(a=a,obj=1) {
  # obj=1 pour plot les valeurs de parametres, obj=2 pour plot les posteriors
    par(mfrow=c(2,2))
    plot(a[[obj]][,1],t="l",ylab="K.pr.Xmin", xlab="Iteration")
    abline(a=0.5,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,5],t="l",ylab="R.pr.Xmin", xlab="Iteration")
    abline(a=0.5,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,2],t="l",ylab="K.pr.Xmax", xlab="Iteration")
    abline(a=10,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,6],t="l",ylab="R.pr.Xmax", xlab="Iteration")
    abline(a=10,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,3],t="l",ylab="K.pr.Xopt", xlab="Iteration")
    abline(a=4,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,7],t="l",ylab="R.pr.Xopt", xlab="Iteration")
    abline(a=4,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,4],t="l",ylab="K.pr.Yopt", xlab="Iteration")
    abline(a=20,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,8],t="l",ylab="R.pr.Yopt", xlab="Iteration")
    abline(a=10,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,9],t="l",ylab="R.tas.Xmin", xlab="Iteration")
    abline(a=270,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,10],t="l",ylab="R.tas.Xmax", xlab="Iteration")
    abline(a=320,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,11],t="l",ylab="R.tas.Xopt", xlab="Iteration")
    abline(a=295,b=0,col="red")
    abline(v=a[[3]],col="blue")

    plot(a[[obj]][,12],t="l",ylab="R.tas.Yopt", xlab="Iteration")
    abline(a=1,b=0,col="red")
    abline(v=a[[3]],col="blue")    
}

plotParallelGibbs <- function(a=a,obj=1) {
    if(obj==1) {
        par(mfrow=c(2,2))
        plot(a[[obj]][1,],t="l",ylab="K.pr.Xmin", xlab="Iteration",col="orange",ylim=c(0,2))
        points(a[[obj+2]][1,],t="l",col="navyblue")

        plot(a[[obj]][5,],t="l",ylab="R.pr.Xmin", xlab="Iteration",col="orange",ylim=c(0,2))
        points(a[[obj+2]][5,],t="l",col="navyblue")        

        plot(a[[obj]][2,],t="l",ylab="K.pr.Xmax", xlab="Iteration",col="orange",ylim=c(8,18))
        points(a[[obj+2]][2,],t="l",col="navyblue")

        plot(a[[obj]][6,],t="l",ylab="R.pr.Xmax", xlab="Iteration",col="orange",ylim=c(8,18))
        points(a[[obj+2]][6,],t="l",col="navyblue")

        plot(a[[obj]][3,],t="l",ylab="K.pr.Xopt", xlab="Iteration",col="orange",ylim=c(2,8))
        points(a[[obj+2]][3,],t="l",col="navyblue")

        plot(a[[obj]][7,],t="l",ylab="R.pr.Xopt", xlab="Iteration",col="orange",ylim=c(2,8))
        points(a[[obj+2]][7,],t="l",col="navyblue")

        plot(a[[obj]][4,],t="l",ylab="K.pr.Yopt", xlab="Iteration",col="orange",ylim=c(15,25))
        points(a[[obj+2]][4,],t="l",col="navyblue")

        plot(a[[obj]][8,],t="l",ylab="R.pr.Yopt", xlab="Iteration",col="orange",ylim=c(8,18))
        points(a[[obj+2]][8,],t="l",col="navyblue")      

        plot(a[[obj]][9,],t="l",ylab="R.tas.Xmin", xlab="Iteration",col="orange",ylim=c(260,280))
        points(a[[obj+2]][9,],t="l",col="navyblue")

        plot(a[[obj]][10,],t="l",ylab="R.tas.Xmax", xlab="Iteration",col="orange",ylim=c(310,330))
        points(a[[obj+2]][10,],t="l",col="navyblue")

        plot(a[[obj]][11,],t="l",ylab="R.tas.Xopt", xlab="Iteration",col="orange",ylim=c(285,305))
        points(a[[obj+2]][11,],t="l",col="navyblue")

        plot(a[[obj]][12,],t="l",ylab="R.tas.Yopt", xlab="Iteration",col="orange",ylim=c(0,2))
        points(a[[obj+2]][12,],t="l",col="navyblue") 


        } else {
            par(mfrow=c(1,1))
            plot(a[[obj]][1,],t="l", ylab="posteriors", xlab="Iteration", col="orange", ylim=c(-1300,-700))
            points(a[[obj+2]][1,], t="l", col="navyblue")
        }
}
