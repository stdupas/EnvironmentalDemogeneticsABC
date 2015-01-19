
### Setting working directories
rm(list=ls())
wd="/media/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # portable
wd="/media/dupas/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # fixe
wd="/home/legs/GraphPOP/" # portable
wd="/home/dupas/GraphPOP/" # fixe
wd="/home/arno/These/GraphPOP" # portable arno
wd="/home/arnaudb/Documents/GraphPOP" # labo arno
setwd(wd)

### Sourcing functions files
source("AskModelsFunctions.R")
source("NicheFunctions.R")
source("DispersionFunctions.R")
source("MutationFunctions.R")
source("CoalescentFunctions.R")

### Sourcing Libraries
library(raster)
library(ape)
library(stringr)
library(lattice)

########### Parameters initialisation  (a priori useless when asking to user...) ########### >>>>>>

###### Environmental data of temperature and precipitations

# Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900)) 
Data2 <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 
# Make raster stack with two layers according to the environmental variables of the dataframe
rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1),"BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)))

###### Genetic parameters :
N=1.5
mutation_rate=1E-4
initial_genetic_value=200

##### Niche Function :

# Concave quadratic skewed distribution parameters
# Shapes of the reaction norms for demographic variables (r and K) :
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
shapesr=c(BIO1="linearPositive",BIO12="conquadraticskewed")
# Parameters of the reaction norm are given by pK and pr
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
#pr = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
pr = matrix(c(100,0.01,NA,NA,NA,NA,NA,NA,NA,NA,300,3000,2500,0,N,N),nrow=8,ncol=2,dimnames=list(c("X0","slope","Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))

# Linear model + constant model parameters
shapesK=c(BIO1="linearPositive",BIO12="linearPositive")
shapesr=c(BIO1="constant",BIO12="constant")
pK = matrix(c(100,0.01,300,0.001),nrow=2,ncol=2,dimnames=list(c("X0","slope"),c("BIO1","BIO12")))
pr = matrix(c(N,N),nrow=1,ncol=2,dimnames=list(c("Y"),c("BIO1","BIO12")))

##### Dispersion Function :

# Shape of the dispersion function :
shapeDisp="fat_tail1"
# Dispersion parameters of the dispersion function
pDisp = c(alpha=1/19,beta=1)

######################### end of parameters initialisation <<<<<<



######################### Using Functions :
##### Get the carrying capacity map :
rasK=rasterStack
values(rasK)= as.matrix(ReactNorm(X=values(rasterStack),p=pK,shapes=shapesK)[,"Y"])
