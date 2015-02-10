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

# load fake GeneticData 
load("GeneticData.RData")

# set the distance method to apply
distanceMethod <- "DeltaMuDistance"

# get the number of simulations (used to initialize reference table)
nbSimul <- 2

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

summaryStatObs <- cbind(ForWLogProb=0,as.data.frame(t(summaryStatObsVect)))



########## Computing summary statistics of simulated data

setwd(paste0(getwd(), "/SimulResults/"))
allFiles <- grep(pattern = "^Genetics_\\d*.txt$", x=list.files(), value = TRUE)

stats <- apply(X = as.array(allFiles), MARGIN = 1, 
               FUN = computeSummaryStats, nbrInd = nbrInd, distanceMethod = "DeltaMuDistance", rotation = rotation)

summaryStatSim = matrix(data = NA, nrow = nbSimul, ncol = ncol(summaryStatObs), dimnames=list(c(), names(summaryStatObs))) 
summaryStatSim[stats[1,],] <- t(stats[-1,])
