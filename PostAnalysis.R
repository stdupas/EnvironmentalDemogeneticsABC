##### POST ANALYSIS

### Setting working directories
rm(list=ls())
wd="/media/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # portable
wd="/media/dupas/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # fixe
wd="/home/legs/GraphPOP/" # portable steph
wd="/home/dupas/GraphPOP/" # fixe
wd="/home/arno/These/GraphPOP" # portable arno
wd="/home/arnaudb/Documents/GraphPOP" # labo arno
setwd(wd)

source("SummaryStats.R")
load("ParamList.RData")
library(abc)

# load fake GeneticData 
GeneticData <- read.table("GeneticData.txt")

# set the distance method to apply
distanceMethod <- "DeltaMuDistance"

# get the number of simulations (used to initialize reference table)
nbSimul <- 10

# number of loci under study:
obsGenetics <- GeneticData[!(colnames(GeneticData)%in%c("x","y","Cell_numbers"))]

# number of individuals under study:
nbrInd <- nrow(obsGenetics)

########## Caracterizing rotation to apply to genetic simulations to get summary statistics

rotation = PCA_rotation(geneticData = obsGenetics, DistanceMethod = distanceMethod)

summaryStatObsMat = as.matrix(do.call(what = distanceMethod, args = list(obsGenetics)))%*%rotation

colN <- paste(rep(colnames(summaryStatObsMat), each = length(rownames(summaryStatObsMat))),
              rownames(summaryStatObsMat),sep=".")

summaryStatObsVect = c(summaryStatObsMat)

names(summaryStatObsVect)<- colN

# First element for package abc
summaryStatObs <- cbind(ForWLogProb=0, as.data.frame(t(summaryStatObsVect)))



########## Computing summary statistics of simulated data

setwd(paste0(getwd(), "/SimulResults/"))
allFiles <- grep(pattern = "^Genetics_\\d*.txt$", x=list.files(), value = TRUE)

stats <- apply(X = as.array(allFiles), MARGIN = 1, 
               FUN = computeSummaryStats, nbrInd = nbrInd, distanceMethod = "DeltaMuDistance", rotation = rotation)

stats <- stats[, order(stats[1,])]

# Second element for package abc
summaryStatSim <- t(stats)[,-1]
colnames(summaryStatSim) <- names(summaryStatObs)

# Third element for package abc :
temp <- as.data.frame(ParamList)
simulParam <- temp[, grep(pattern = "Values", x = names(temp))]

########### ABC package... Let's go giiiirls !

abc(target = summaryStatObs, param = simulParam, sumstat = summaryStatSim, tol = 0.5, method = "loclinear" )



