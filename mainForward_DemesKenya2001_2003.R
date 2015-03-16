### Setting working directories
rm(list=ls())
wd="/media/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # portable
wd="/media/dupas/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # fixe
wd="/home/legs/GraphPOP/" # portable steph
wd="/home/dupas/GraphPOP/" # fixe
wd="/home/arno/These/GraphPOP" # portable arno
wd="/home/arnaudb/Documents/GraphPOP" # labo arno
setwd(wd)

### Sourcing functions files
source("AskModelsFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
#source("MutationFunctions.R")
#source("CoalescentFunctions.R")
source("PriorFunctions.R")
#source("MarkovProcess.R")
#source("GeneticDataSimulation.R")
source("forwardModels.R")
source("readNetCDF.R")

### Sourcing Libraries
library(raster)
library(ape)
library(stringr)
library(lattice)
library(parallel)
library(lubridate) # of little use, should be removed (function "days()")
library(RNetCDF) # to read NetCDF data (climate spatial time series)

########### Parameters initialisation  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

###### Environmental data of temperature and precipitations
# Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900)) 
#Data2 <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 
EnvDataRasterStack = nc2EnvDataAndRasterStack(ncDirectory=paste(wd,"ForwardSimulData/",sep=""),aggregationParam=10)
#EnvDataRasterStack = nc2EnvDataAndRasterStack(ncDirectory=paste(wd,"ForwardSimulData/",sep=""),aggregationParam=4)
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

minDates = min(release$birthDate)
maxDates = min(max(recovery$birthDate),max(as.Date(colnames(EnvData))))
Dates <- as.Date(as.Date(minDates):as.Date(maxDates),origin="1970/01/01")

recovery <- read.table("ForwardSimulData/Stemborer_Kenya2001_2005.csv")
recovery <- recovery[,c("Long_dec","Lat_dec","Diss_Date","B._fusca")]
colnames(recovery) <- c("x","y","birthDate","size")
recovery$birthDate <- as.Date(as.character(recovery$birthDate),format="%d/%m/%y")
recovery$demeNb <- cellFromXY(object = rasterStack, xy = recovery[, c("x", "y")])
recovery <- aggregation(recovery,BY=c("birthDate","demeNb"),methodes=c("Mean","Mean","Name","Sum","Name"))
recovery <- recovery[which(as.Date(recovery$birthDate)<=as.Date("2003/12/31")),]

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

distMat <- distanceMatrixFromRaster(object = rasterStack)*111.32
demeSizes <- array(0,dim=c(nrow(EnvData),length(Dates)),dimnames = list(1:nrow(EnvData),as.character(Dates)))
demeSizes[,as.character(birthDates)] <- 1
EnvData <- EnvData[,colnames(demeSizes),]

likelihoodShort()

nlm("likelihoodShort",p=c(dispersionRate = .025,dispersionDistance=100,
                          K.pr.X0=-5,K.pr.Xopt=267.1267,K.pr.Yopt=10,K.generationTime=25,K.generationTimeSD=3,
                          R.pr.X0=-5,R.pr.Xopt=267.1267,R.pr.Yopt=10,R.generationTime=25,R.generationTimeSD=3))
