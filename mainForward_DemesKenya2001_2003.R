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
library(bbmle) # Library fr the function of minimization mle2

########### Parameters initialisation  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

###### Environmental data of temperature and precipitations
# Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900)) 
#Data2 <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 

#EnvDataRasterStack = nc2EnvDataAndRasterStack(ncDirectory="/Users/Stagiaire/ForwardSimulData/",aggregationParam=10)
EnvDataRasterStack = readRDS("/Users/Stagiaire/ForwardSimulData/ObjectEnvdataRasterStack")

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

########################## VERIFIER FORMAT FICHIER CSV #########################################
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

recovery2 = buildDataSet()

for ( i in 1:dim(EnvData2)[1]){
    for (j in 1:dim(EnvData2)[2]){
        for (k in 1:dim(EnvData2)[3]){
            if(is.na(EnvData2[i,j,k]) == TRUE){
                EnvData2[i,j,k] = 0
            }
        }
    }
}



likelihoodShort()


test = nlm(f = likelihoodShortTest,p=c(dispersionRate = .025,dispersionDistance=1,
                                       K.pr.X0=0,K.pr.Xopt=38.40947,K.pr.Yopt=11.53846,
                                       R.pr.X0=0,R.pr.Xopt=38.40947,R.pr.Yopt=1,
                                       generationTime=25,generationTimeSD=3, 
                                       dvlpTime=5, dvlpTimeSD=1))

testShort = nlm(f = likelihoodShortTest,p=c(K.pr.X0=5,K.pr.Xopt=30,K.pr.Yopt=15,
                                            R.pr.X0=5,R.pr.Xopt=30,R.pr.Yopt=3),
                print.level=2, ndigit=3)
