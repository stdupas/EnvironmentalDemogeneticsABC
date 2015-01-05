####
#### STATISTIC MODEL : INFERENCE OF ECOLOGICAL MODEL FROM GENETIC DATA
####

##########################################################################
############## Set your working directory and files to load ##############
##########################################################################

# matrice : colonne
#raster : ligne


#############################
#### Library to download ####
#############################

library(rgdal)
library(raster)
library(MASS)
library(Geneland)
library(ape)
library(stringr)
library(lattice)
#library(RNetCDF)

####################################################
##### BACKWARD MODEL FUNCTIONS AND EXECUTIONS ######
#####  SIMULATION OF PREDICTED GENETIC DATA   ######
####################################################


# Function to get spatial resolution in km
degree2km = function(rasterStack){
  x_origin = ((xmin(rasterStack)+xmax(rasterStack))/2) #longitude origin
  y_origin = ((ymin(rasterStack)+ymax(rasterStack))/2) #latitude origin
  x_destination = (x_origin + xres(rasterStack)) #longitude of destination point
  y_destination = (y_origin + yres(rasterStack)) #latitude of destination point
  
  dist_degree <- acos(sin(x_origin)*sin(x_destination)+cos(x_origin)*cos(x_destination)*cos(y_origin-y_destination))
  dist_km = dist_degree * 111.32
  dist_km
}

# Aggregate_and_adjust_raster_to_data change resolution and extent of environmental stacked layers
# according to data geographic range and extension zone outside geographic range of data
# ARGUMENTS:
# Envir_raster_stack = raster file 
# release = release points file (columns "X" and "Y" as longitude nd latitude)
# recovery = recovery points file (columns "X" and "Y" as longitude nd latitude)
Aggregate_and_adjust_raster_to_data <- function(Envir_raster_stack,release,recovery,extend_band_size,aggregate_index)
{
  samples <- SpatialPoints(rbind(na.omit(release[,c("X","Y")]),na.omit(recovery[,c("X","Y")])))
  if (aggregate_index > 1) {Envir_raster_stack <- aggregate(crop(Envir_raster_stack,extent(samples)+extend_band_size), fact=aggregate_index, fun=mean, expand=TRUE, na.rm=TRUE)} else {
    Envir_raster_stack <- crop(Envir_raster_stack,extent(samples)+extend_band_size)
  }
  Envir_raster_stack
}

################################################
# ENVIRONMENTAL NICHE FUNCTIONS
# The following functions are used to transform environmental variables to 
# growth rate and carrying capacity
#
# We have to define specific reaction norms for these population traits
# 1/ log.additive - the first R and K reaction norm model for several environmental variables 
# is the log additive
# R and K depend on the product of reaction norms for each variable
# for instance assume insufficient precipitationq reduces number of offspring by 0.8
# and inadequate temperature by 0.6, overall the number of offspring will be
# reduced by 0.6 *0.8 = 0.48. In this model we assume that the two factors acts independently 
# on log values of each variable. on log scale reaction norms are additive
# 2/ 
# conquadraticiopt is an asymetric conquadratic density response fonction 
#Xmax=p[2];p["Xmin"]=p[1];Ymax=p["Ymax"];p["Xopt"]=p["Xopt"]


# conquadraticskewed
# asymetric concave conquadratic function 
# X : matrix or data frame providing the values of independent variable to calculate reaction norm
# p : matrix parameter values for the reaction norm
# line names of p : c("Xmin","Xmax","Xopt","Ymax"), 
# column names of p : names of the independent variables of X used for the reaction norm calculation
# [Xmin, Xmax] is the enveloppe, 
# Xopt is the value that maximises the function
# Yopt is the maximum value of the function

conquadraticskewed <- function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))))
{
  Yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]
  Xopt = p[rep("Xopt",dim(X)[1]),colnames(X)]
  Xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]
  Xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
  alpha <- -log(2)/log((Xopt-Xmin)/(Xmax-Xmin))
  Xprime<- ((X-Xmin)/(Xmax-Xmin))^alpha*(Xmax-Xmin)+Xmin
  y <- (Yopt-4*Yopt/((Xmax-Xmin)^2)*(Xprime-(Xmin+Xmax)/2)^2)*((X>=Xmin)&(X<=Xmax))
  y[X<Xmin] <-0
  y
}

conquadraticskewedsq <- function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))))
{
  Yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]
  Xopt = p[rep("Xopt",dim(X)[1]),colnames(X)]
  Xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]
  Xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
  alpha <- -log(2)/log((Xopt-Xmin)/(Xmax-Xmin))
  Xprime<- ((X-Xmin)/(Xmax-Xmin))^alpha*(Xmax-Xmin)+Xmin
  y <- (Yopt-4*Yopt/((Xmax-Xmin)^2)*(Xprime-(Xmin+Xmax)/2)^2)*((X>=Xmin)&(X<=Xmax))
  y+y*(Yopt-y)
#  y[X<Xmin] <-0
}

# conquadratic = concave quadratic function 
# X : matrix or data frame providing the values of independent variable to calculate reaction norm
# p : matrix parameter values for the reaction norm
# line names of p : c("Xmin","Xmax","Xopt","Ymax"), 
# column names of p : names of the independent variables of X used for the reaction norm calculation
# [Xmin, Xmax] is the enveloppe, 
# Xopt, is not considered and is given the value (Xmin+Xmax)/2
# Yopt is the maximum value of the function

conquadratic <- function(X,p)
{
  xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
  xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]  
  yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]  
  (yopt-(4*yopt/(xmax-xmin)^2)*(X-(xmin+xmax)/2)^2)*((X>xmin)&(X<xmax))
}

# conquadraticqsq = concave squared quadratic function 
# X : matrix or data frame providing the values of independent variable to calculate reaction norm
# p : matrix parameter values for the reaction norm
# line names of p : c("Xmin","Xmax","Xopt","Ymax"), 
# column names of p : names of the independent variables of X used for the reaction norm calculation
# [Xmin, Xmax] is the enveloppe, 
# Xopt, is not considered and is given the value (Xmin+Xmax)/2
# Yopt is the maximum value of the function

conquadraticsq <- function(X,p)
{
  xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
  xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]  
  yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]  
  res = (yopt-(4*yopt/(xmax-xmin)^2)*(X-(xmin+xmax)/2)^2)*((X>xmin)&(X<xmax))
  res + res * (1-res)
}

# envelope = simple envelope function 
# X : matrix or data frame providing the values of independent variable to calculate reaction norm
# p : matrix parameter values for the reaction norm
# line names of p used for caclculation : c("Xmin","Xmax","Yopt"), 
# column names of p : names of the independent variables of X used for the reaction norm calculation
# [Xmin, Xmax] is the enveloppe, 
# Yopt is the value of the function in the [Xmin, Xmax] interval

enveloppe <- function(X,p)
{
  p[rep("Yopt",dim(X)[1]),colnames(X)]*((X>p[rep("Xmin",dim(X)[1]),])&(X<p[rep("Xmax",dim(X)[1]),colnames(X)]))
}

# envelope = linear response within an envelope
# X : matrix or data frame providing the values of independent variable to calculate reaction norm
# p : matrix parameter values for the reaction norm
# line names of p used for caclculation : c("Xmin","Xmax","Yxmin","Yxmax"), 
# column names of p : names of the independent variables of X used for the reaction norm calculation
# [Xmin, Xmax] is the enveloppe, 
# Yxmin and Yxmax are the values at Xmin and Xmax

envelinear <- function(X,p,log=FALSE)
{
  Yxmin = p[rep("Yxmax",dim(X)[1]),colnames(X)]
  Yxmax = p[rep("Yxmin",dim(X)[1]),colnames(X)]
  Xmin=  p[rep("Xmin",dim(X)[1]),colnames(X)]
  Xmax=  p[rep("Xmax",dim(X)[1]),colnames(X)]
  a = (Yxmin - Yxmax) / (Xmin - Xmax)
  b = Yxmin - Xmin * a
  if (log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
}
  
envelin0 <- function(X,p,log=FALSE)
{
  Yxmin = p[rep("Yxmax",dim(X)[1]),colnames(X)]
  Yxmax = p[rep("Yxmin",dim(X)[1]),colnames(X)]
  Xmin=  p[rep("Xmin",dim(X)[1]),colnames(X)]
  Xmax=  p[rep("Xmax",dim(X)[1]),colnames(X)]
  a = (Yxmin - Yxmax) / (Xmin - Xmax)
  b = Yxmin - Xmin * a
  if (log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
}

linear <- function(X,p)
{
  Yx1 = p[rep("Yx1",dim(X)[1]),colnames(X)]
  Yx0 = p[rep("Yx0",dim(X)[1]),colnames(X)]
  a = (Yx1 - Yx0)
  b = Yx0
  if (log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
}

#inc=(X-xmin)/(Xmax-xmin)*p["Ymax",]*(X>xmin)*(X<xmax),
#lineardec=(xmax-X)/(Xmax-xmin)*p["Ymax",]*(X>xmin)*(X<xmax),


#bi <- mulconquadraticskewed(Data2,p=matrix(c(100,400,250,1,300,2000,1000,1),nrow=4,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Ymax"),c("BIO1","BIO12"))))
#conquadraticskewed[]
#(p[4]-4*p[4]/((p[2]-p[1])^2)*(( ((X-p[1])/(p[2]-p[1]))^-log(2)/log((p["Xopt"]-p[1])/(p[2]-p[1]))*(p[2]-p[1])+p[1])-(p[1]+p[2])/2)^2)*(X>=p[1])*(X<=p[2])

# ReactNorm computes reaction norm value and geometric mean of the reaction 
# norms for each variable in a data frame
# p are the parameter values of the reaction norm for each variable
# shape is the shape of the reaction norm

ReactNorm <- function(X,p,shapes)
  #p=c(p["Xmin"]=10,p["Xmax"]=20,p["Xopt"]=18,p["Ymax"]=0.1)
  ## shapeDisps in c("enveloppe","envelin","envloglin","loG","conquadratic","conquadraticskewed","conquadraticsq","conquadraticskewedsq")
{
  Y=X
  if (!all(colnames(p)%in%names(shapes))) {stop ("variable names do not correspond between parameters 'p' ans 'shape'")}
  for (shape in as.character(levels(as.factor(shapes))))
  {
    variables = colnames(p)[which(shapes==shape)]
    Y[,variables]=switch(shape,
           enveloppe=enveloppe(subset(X,select=variables),p),
           envelin=envelinear(subset(X,select=variables),p),
           envloglin=envelinear(subset(X,select=variables),p,log=TRUE),
           loG = log(subset(X,select=variables)),
           linear = linear(subset(X,select=variables),p),
           conquadratic=conquadratic(subset(X,select=variables),p),
           conquadraticskewed=conquadraticskewed(subset(X,select=variables),p),
           conquadraticsq=conquadraticsq(subset(X,select=variables),p),
           conquadraticskewedsq=conquadraticskewedsq(subset(X,select=variables),p)
    )
  }
  Y=cbind(Y,Y=apply(Y, 1, prod)^(1/dim(p)[2])) # geometric mean
  Y
}


# K_Function: Gets effective carrying capacity from environmental variables 
# p : parameters 
# rasterStack = environmental variables 
# shapes = vector of shapes used for the niche function of each environemental variable 
# (among enveloppe, envelin, envloglin, conquadratic, conquadraticskewed) 
K_Function <- function(rasterStack, p, shapes){
  #K=matrix(combineReactNorms(values(rasterStack),p,shapes),byrow=TRUE,nrow=length(rasterStack),ncol=legnth(rasterStack))
  ReactNorm(values(rasterStack),p,shapes)
}

# R_Function: Gets effective growth rate from environmental variables 
# alpha and beta are fixed by estimation
# rasterStack = environmental variables 
R_Function <- function(rasterStack, alpha, beta){
  if(nlayers(rasterStack)>1){
    R = exp(as.matrix(alpha+sum(beta*rasterStack)))# utilisation d'un modele lineaire generalise
  }
  else{ R = exp(as.matrix(alpha+beta*rasterStack)) }
  R = t(R) # transpose to get niche predicted values that fits to matrix organisation (by columns) and not in raster organisation (by rows)
  R = t(matrix(R,nrow=length(R),ncol=length(R))) # Get population size by columns
  R[is.na(R)]<-0 # replace NA by 0
  R
}

Show_Niche <- function(Data,p,shapes=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")) # non terminé
{
  pairs = NULL;i=0
  while (i < dim(p)[2])
  {
    i=i+1;j=i+1
    first= colnames(p)[i]
    while (j <= dim(p)[2])
    {
      second = colnames(p)[j]
      j=j+1
      pairs=rbind(pairs,c(first,second))
    }
  }
  pairs
  dev.off()
  par(mfrow=c(1,1))
  Data= as.data.frame(Data)
  for (i in 1:dim(pairs)[1])
  {
    print(i)
    form = as.formula(paste("z~",paste(pairs[i,],collapse="*"),sep=""))
    Data[,"z"]=ReactNorm(Data,p,shapes)[,"Y"]
    wireframe(form,data=Data,scales=list(arrows=FALSE)) # requires library lattice
  }
}

# populationSize: uses K_Function to obtain the populationsize landscape raster
#
#
populationSize <- function(donneesEnvironmentObs, p, shapes)
{
  populationSize <- donneesEnvironmentObs
  values(populationSize) <- ReactNorm(valules(donneesEnvironmentObs), p, shapes)[1,]
  populationSize
}
  
# (optional) distanceMatrix return distance between all cells of raster
distanceMatrix <- function(rasterStack){
  #get x and y coordinates for each cell of raster object put in parameters
  coords = xyFromCell(rasterStack, 1:length(values(rasterStack[[1]])), spatial=FALSE)
  distance = as.matrix(dist(coords)) # distance matrix of coordinates
  return(distance)
}

# migrationMatrix return matrix of migration rate for each cell
# From the raster Stack and the dispersal parameter pDisp (estimated),
# this function calculates distance between all cells of raster
# To this matrix of distance is applied a shapeDisp of migration.
# And this function allows to choose the shapeDisp of study:
# shapeDisp:  "gaussian" a simple normal density distribution, 
#           "exponential" density distribution,
#           "fat_tail1", "fat_tail2": ref :Chapman et all, Journal of Animal Ecology (2007) 76 , 36– 44
#           "island" probability 1-m to stay, else homogen dispersion,
#           "contiguous" near dispersal.
migrationMatrix <- function(rasterStack,shapeDisp, pDisp){
  coords = xyFromCell(rasterStack, 1:length(values(rasterStack[[1]])), spatial=FALSE)
  distanceMatrix = as.matrix(dist(coords)) 
  migration = apply(distanceMatrix, c(1,2), 
                    function(x)(switch(shapeDisp,
                                       # 1: alphaDisp   2: betaDisp ; note: esperance = 1/alphaDisp
                                       fat_tail1 = 1/(1+x^pDisp[2]/pDisp[1]), # Molainen et al 2004
                                       # 1: sigmaDisp
                                       gaussian = (dnorm(x, mean = 0, sd = pDisp[1], log = FALSE)),
                                       # 1: sigma Disp
                                       exponential = (dexp(x, rate = 1/pDisp[1], log = FALSE)),
                                       # 1: sigmaDisp
                                       contiguous = (x==0)*(1-pDisp[1])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp[1]/2),
                                       #island = (x==0)*(1-pDisp[1])+(x>0)*(pDisp[1]/(nCell-1)),
                                       island = (x==0)*(1-pDisp[1])+(x>0)*(pDisp[1]),
                                       #1: sigmaDisp    2: gammaDisp
                                       fat_tail2 = x^pDisp[2]*exp(-2*x/(pDisp[1]^0.5))
                    )))
  return(migration)
}

# transitionMatrix obtained with an isotropic migration hypothesis for a backward model
transitionMatrixBackward <- function(r,K, migration){
  if ((length(r)==1)&(length(K)==1)){transition = r * K * t(migration)}
  if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
  if ((length(r)==1)&(length(K)>1)){transition = r * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
  if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
  if ((length(r)>1)&(length(K)>1)) {
  transition = t(matrix(r,nrow=length(r),ncol=length(r))) * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)
  transition = transition / rowSums(transition)  # Standardisation
  transition = transition * as.numeric(transition>1e-6) # removal of values below 1e-6
  transition = transition / rowSums(transition)  # Standardisation again
  }
  transition
}

# transitionMatrix obtained with an isotropic migration hypothesis for a formard model
transitionMatrixForward <- function(r,K, migration, meth="non_overlap"){
  rs = matrix(r,nrow=length(r),ncol=length(r))
  Ku = t(matrix(K,nrow=length(K),ncol=length(K)))
  leave = migration*(1+rs)*t(Ku); leave = leave - diag(leave)
  switch (meth,
  non_overlap = migration * rs * Ku / colSums(rs * t(Ku) * migration),
  overlap = migration * (1+rs) * Ku / (colSums((1+rs) * t(Ku) * migration - t(leave)))
  )
}  
#  {
#    transition = migration * Npop # number of individuals going from cell i to cell j
#    transition = transition / t(matrix(rowSums(transition),nrow=nCell,ncol=nCell))  # Standardisation by total number of individuals reaching cell j
#    transition
#  }
#}

#prior

# simulation forward of population sizes across time
forward_simul_landpopsize <- function(N0,p, migration)
{
  
}
# laplaceMatrix returns Laplacian matrix from transition matrix
laplaceMatrix <- function(transitionMatrix){
  matrixD = matrix(0,nrow = dim(transitionMatrix), ncol = dim(transitionMatrix))
  diag(matrixD) = 1 # diagonal equals to 1
  laplacianMatrix = matrixD - transitionMatrix
  laplacianMatrix[is.na(laplacianMatrix)]<-0 # replace NA by 0
  laplacianMatrix
  
}

# Calcul of resistance between two points of the graph
# with the Moore-Penrose generalized inverser matrix.
# ref : Bapat et all, A Simple Method for Computing Resistance Distance (2003)
# ref : Courrieu, Fast Computation of Moore-Penrose Inverse Matrices (2005)
resistDist <- function(laplacianMatrix){
  inverseMP = ginv(laplacianMatrix) # generalized inverse matrix  (Moore Penrose)
  diag = diag(inverseMP) # get diagonal of the inverse matrix
  mii = matrix(diag, nrow =dim(inverseMP), ncol = dim(inverseMP))
  mjj = t(mii)
  mij = inverseMP
  mji = t(mij)
  commute_time = mii + mjj - mij - mji
  commute_time
}


# Calcul of genetic distance from resistance 
geneticDist <- function(commute_time, popSize){
  #genetic_dist = commute_time / (8* popSize)
  genetic_dist = commute_time / (8* (sum(popSize)/(dim(popSize)[1]*dim(popSize)[2])))
  genetic_dist
}

# MAIN EXECUTIONS DES FONCTIONS PERMETTANT D'OBTENIR IN FINE LES DONNEES GENETIQUES PREDITES
#rasterCrop = Aggregate_and_adjust_raster_to_data(raster(paste(wd,envdir,envfiles,sep="")),release=read.table(paste(wd,genetfile,sep="")), recovery=read.table(paste(wd,genetfile,sep="")), extend_band_size=1, aggregate_index=aggregate_factor)
#plot(rasterCrop)
  

###################################################
##### FORWARD MODEL FUNCTIONS AND EXECUTIONS ######
#####  SIMULATION OF OBSERVED GENETIC DATA   ######
###################################################


# This function returns a data frame that contains genetic data:
# coords X and Y of each cell - number of cells - locus with 2 alleles
# ARGUMENTS: 
# RasPopSizes = raster of pop sizes
# dataSize = data that contains number of individuals in each cell
# nb_locus = number of locus wished
# initial_locus_value = genetic information for the first generation (parents)
# option =  sample/full_1col/2col_haploid/diploid
#           sample : nind individuals are sampled from the full population
#           full : number of individuals per cell equals carrying capacity
#           1col/2col : for diploids, one columns or 2 columns representation
#           haploid or diploid : haploid or diploid 
# nind = number of individuals sampled from the full carrying capacity population

CreateGenetArray <- function(rasK, nb_locus, initial_locus_value,Option="sample_1_col_diploid",nind=4)
{
  #Get coords for each cell
  coords = xyFromCell(rasK, 1:length(values(rasK[[1]])), spatial=FALSE)
  repet = switch(Option,
                 sample_1col_diploid = sort(rep(sample(rep(1:ncell(rasK),round(values(rasK))),nind),2)), 
                 sample_2col_diploid = sample(rep(1:length(rasK),round(values(rasK))),nind), 
                 sample_haploid = sample(rep(1:ncell(rasK),round(values(rasK))),nind),
                 full_1col_diploid = rep(1:length(rasK),round(values(rasK))*2), 
                 full_2col_diploid = rep(1:length(rasK),round(values(rasK))), 
                 full_haploid = rep(1:length(rasK),round(values(rasK))),                 
  )
  geneticData <- as.data.frame(coords[repet,]) ; geneticData[,"Cell_numbers"] <- repet
  genes = switch(Option,
                 sample_1col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                 sample_2col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                 sample_haploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                 full_1col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                 full_2col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                 full_haploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),                 
  )
  colnames(genes) = switch(Option,
                           sample_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                           sample_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                           sample_haploid = paste("Locus",1:nb_locus, sep=""),
                           full_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                           full_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                           full_haploid = paste("Locus",1:nb_locus, sep=""),                 
  )
  geneticData <- cbind(geneticData,genes) # add locus to geneticData 
  geneticData
}



#Plot genetic data in environmental data observed 
#Genetic data is turn into a Spatial Pixel Data Frame 
#Mettre des couleurs en fonction du nombre d'individu
plotGeneticData = function(geneticData, EnvironmentalDataObserved){
  colnames(geneticData)[1:2] <- c("x","y")
  geneticData = SpatialPixelsDataFrame(points = geneticData[,c("x","y")], data = geneticData[,])
  plot(EnvironmentalDataObserved[[1]])
  plot(geneticData, add = T)
}

# Function that computes population size distribution moments in a grid from one generation to the other
# N population sizes of the parent genration
# r growth rates
# K carrying capacities
# d death rates
# migration transition matrix between cells of the grid
# ptG : parameters of generaiton time model

gridRepnDispFunction <- function(dynamics,r,K,d=.9,ptG, migration,overlapping=TRUE)
{
  # values(dynamics)[,dim(values(dynamics))[2]] is value at previoud day
  # d is mortality
  Nt = values(dynamics)[,dim(values(dynamics))[2]]*(1-d) + r*N*(K-N/K)
  esperance[K==0] <- 0  
}

# Function that combine reproduction, dispersion and mutation for a given genetic data
repnDispMutFunction <- function(geneticData, dimGeneticData, mutationRate, transitionmatrice){
  # Calcul for reproduction and dispersion 
  # random choice of individuals and at the same time of their target cells in the transition matrix
  ncell_transition <- dimGeneticData[1]*dim(transitionmatrice)[1]
  transition_celnu <- sample(ncell_transition,dimGeneticData[1],replace=FALSE,transitionmatrice[geneticData[,"Cell_numbers"],])
  transition_col <- ceiling(transition_celnu/dimGeneticData[1])
  transition_line <- transition_celnu%%(dimGeneticData[1]);transition_line[transition_line==0]<-dimGeneticData[1]
  cell_numbers_sampled <- geneticData[transition_line,"Cell_numbers"]
  geneticData <- geneticData[transition_line,]
  geneticData[,"Cell_numbers"] <- transition_col
  locusCols = grep("Locus", colnames(geneticData))
  step = 2
  mu = mutationRate # mutation rate
  liability = runif(prod(dimGeneticData), 0, 1) # mutation liability
  liability = as.data.frame(matrix(liability, ncol = length(locusCols), nrow = dimGeneticData[1]))
  geneticData[,locusCols] = geneticData[,locusCols] + ((liability<mu/2)*step - (liability>(1-mu/2))*step) 
  #print(c("mutrat",(sum(liability<mu/2)+sum(liability>(1-mu/2)))/(length(grep("Locus", colnames(geneticData)))*dimGeneticData[1])))
  geneticData
}

# Function that combine reproduction, dispersion and mutation for a given genetic data
repnDispMutFunction <- function(geneticData, dimGeneticData, mutationRate, transitionmatrice){
  # Calcul for reproduction and dispersion 
  # loop for individuals
  locusCols = grep("Locus", colnames(geneticData))
  for (individual in 1:dimGeneticData[1])
  { # we choose where the parent come from in the individual cell line probabilities of the backward transition matrice
    mothercell = sample(nCell, 1,,transitionmatrice[geneticData[individual,"Cell_numbers"],])
    # we chose the parent among the individuals in this cell
    geneticline = sample(which(geneticData[,"Cell_numbers"]==mothercell),1)
    # we atribute the individual the genetic data of its mother
    geneticData[individual,locusCols] = geneticData[geneticline,locusCols]
  }
  step = 2
  mu = mutationRate # mutation rate
  liability = runif(prod(dimGeneticData), 0, 1) # mutation liability
  liability = as.data.frame(matrix(liability, ncol = length(locusCols), nrow = dimGeneticData[1]))
  geneticData[,locusCols] = geneticData[,locusCols] + ((liability<mu/2)*step - (liability>(1-mu/2))*step) 
  geneticData
}

#Function that calculate probability of identity of genes intra individual (at individuals level)
Qwithin_pair <- function(geneticData){
  matrix_pair = geneticData[,grep(".2", colnames(geneticData), fixed = T)]
  matrix_impair = geneticData[,grep(".1", colnames(geneticData), fixed = T)]
  Qw <- (matrix_pair == matrix_impair)
  Qw = rowMeans(Qw) # vector of probability of Qw for each individual
  Qw = matrix(Qw, ncol = dim(geneticData)[1], nrow = dim(geneticData)[1])
  Qw = (Qw+t(Qw))/2
  Qw
}

#Function that calculate probability of identity of genes intra individual (at population level)
Qwithin_pop <- function(geneticData){
  matrix_pair = geneticData[,grep(".2", colnames(geneticData), fixed = T)]
  matrix_impair = geneticData[,grep(".1", colnames(geneticData), fixed = T)]
  Qw <- (matrix_pair == matrix_impair)
  Qw = mean(Qw*1)
  Qw
}


# Fonction that calculates probability of identity of genes inter individual (between two individuals)
Qbetween <- function(geneticData, dimGeneticData){
  Qb = matrix(ncol = dimGeneticData[1], nrow = dimGeneticData[1]) #initialization of Qb as a matrix
  # A = genetic data with loci only
  A=as.matrix(geneticData[,grep("Locus",colnames(geneticData),fixed=T)])
  # On construit un tableau à plusieurs étages: il y a autant d'étage qu'il y a de locus ( les alleles passent en étage)
  A3 = aperm(array(A,dim=c(dim(A)[1],dim(A)[2],dim(A)[1])), c(1,3,2)) # permutation des colonnes et des etages
  B3 = aperm(A3, c(2,1,3)) # transposee de A3
  moy1 = colMeans(aperm(A3 == B3), dims = 1, na.rm = T)
  
  #Permutation des colonnes deux à deux des alleles de A pour calculer l'autre cas possible d'identite des alleles
  l= 1:dim(A)[2]
  Aprime= A[,c(matrix(c(l[2*floor(1:(length(l)/2))],l[2*floor(1:(length(l)/2))-1]), nrow= 2, ncol = length(l)/2, byrow = T))] # Permutation des colonnes
  # permutation et creation des etages pour Aprime ( on ne change qu'une seule des matrice (A3 / B3) et l'autre est inchangée pour comparer les alleles"complementaires"
  Aprime3 = aperm(array(Aprime,dim=c(dim(A)[1],dim(A)[2],dim(A)[1])), c(1,3,2))
  moy2 = colMeans(aperm(Aprime3 == B3), dims = 1, na.rm = T) # calcul moy dist pour les individus avec loci permutés Aprime3 et la transposée B3
  #Mean of distance between individuals: Qbetween
  Qb =(moy1 + moy2)/2
  Qb
}


# TEST: pour un nombre de generation donné, on teste la stabilité du a value
#Function test of stabilisation for a value
test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice){
## ref: Rousset et al. J Evol Biol,13  (2000) 58-62.
  vecteur_a_value <-c(0)
  for(i in 1: nb_generations){
    print(i)
    geneticData <- repnDispMutFunction(geneticData, dimGeneticData, mutationRate, transitionmatrice) 
    matrixQb = (1-Qwithin_pair(geneticData)+2*(Qwithin_pair(geneticData)-Qbetween(geneticData,dimGeneticData)))
    matrixQw = 2*(1-Qwithin_pop(geneticData))
    vecteur_a_value[i] = matrixQb/matrixQw-1/2
    vecteur_a_value[is.na(vecteur_a_value)] <-0
    if((i>90) && (i%%30 == 0)){
      if(var(vecteur_a_value[(i-30):i])> var(vecteur_a_value[(i-60):(i-30)]) 
         && var(vecteur_a_value[(i-30):i])> var(vecteur_a_value[(i-90):(i-60)])){
        return(list(geneticData, (matrixQb/matrixQw-1/2)))
        break
      }
    }
  }
}

test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice){
  ## ref: Rousset et al. J Evol Biol,13  (2000) 58-62.
  vecteur_a_value <-c(0)
  for(i in 1: nb_generations){
    print(i)
    geneticData <- repnDispMutFunction(geneticData, dimGeneticData, mutationRate, transitionmatrice) 
    Genotypes = geneticData[,grep("Locus", colnames(geneticData), fixed = T)]
    popmbrship=geneticData[,"Cell_numbers"]
    Fst = Fstat(Genotypes,nCell,popmbrship,ploidy=2)
    Fstlinear[i] <- Fst/(1-Fst)
    if((i>90) && (i%%30 == 0)){
      if(var(Fst[(i-30):i])> var(Fst[(i-60):(i-30)]) 
         && var(Fst[(i-30):i])> var(Fst[(i-90):(i-60)])){
        return(list(geneticData, Fst))
               break
      }
    }
  }
}



fstat = function(geneticData){
  Genotypes = geneticData[,grep("Locus", colnames(geneticData), fixed = T)]
  form <- as.formulae
  Pops = geneticData[,"Cell_numbers"]
  MeanPop = t(matrix((colSums(Genotypes)/dimGeneticData[1]),ncol=dimGeneticData[1],nrow=dimGeneticData[2]))
  VarInd = matrix(Genotypes^2 - MeanTot^2,ncol=dimGeneticData[2])
  VarInterPop = var(MeanPop)
  VarIntraPop = colSums(VarInd)/dimGeneticData[1]
  VarTot = VarInd
}

simul_commute <- function(cells=c(1,2),transitionmatrice)
{
  tmpcell <- cells[1];t=1
  while (tmpcell != cells[2])
  {
    tmpcell =  sample(dim(transitionmatrice)[2],size=1,prob=c(transitionmatrice[tmpcell,]))
    t=t+1
  }
  hit=TRUE
  while (tmpcell != cells[1])
  {
    tmpcell =  sample(dim(transitionmatrice)[2],size=1,prob=c(transitionmatrice[tmpcell,]))
    t=t+1
  }
  commute=TRUE
  t
}

simul_coocur <- function(cells=c(1,2),transitionmatrice)
{
  tmpcell1 <- cells[1];tmpcell2 <- cells[2];t=1
  while (tmpcell1 != tmpcell2)
  {
    tmpcell1 =  sample(dim(transitionmatrice)[2],size=1,prob=c(transitionmatrice[tmpcell1,]))
    tmpcell2 =  sample(dim(transitionmatrice)[2],size=1,prob=c(transitionmatrice[tmpcell2,]))
    t=t+1
  }
t
}

# simul_coalescent
# function that simulates a coalescent in a lansdcape characterized by an environmental variable raster
# stack in aspecies with given niche function relating the carrying capacity and growth rate 
# with the environemental variable 
# pK : parameters of K / environmental variables
# pr : parameters of r / environmental variables
# shapesK : shapes of niche funciotn (reaction norm)
# pDisp : parameters of dispersion 
# rasterStack : environmental variables raster stack
# genetic data table : with coordinates
# aggre_gener : number of generations to aggregate in the simulation steps

simul_coalescent <- function(geneticData,rasterStack,pK,pr,shapesK,shapesr,shapeDisp,pDisp,mutation_rate=1E-4,aggre_gener=1E3)
{
  prob_forward=NA
  K = ReactNorm(values(rasterStack),pK,shapesK)[,"Y"]
  r = ReactNorm(values(rasterStack),pr,shapesr)[,"Y"]
  migrationM <- migrationMatrix(rasterStack,shapeDisp, pDisp)
  transitionmatrice = transitionMatrixBackward(r, K, migration= migrationM);
  transition_forward = transitionMatrixForward(r, K, migration= migrationM)
  N <- round(K)
  coalescent = list() # list containing all the times and genes conserved by coalescent events
  # when 2 genes or more coalesce, only the the genes tagged by has the smallest number remains
  nodes = as.numeric(rownames(geneticData));names(nodes)=as.character(nodes)# names of the tip nodes that will coalesce
  cell_number_of_nodes <- geneticData[,"Cell_numbers"] # where were the genes sampled in the landscape
  names(cell_number_of_nodes) <- nodes
  parent_cell_number_of_nodes <- cell_number_of_nodes # where the previous generation genes were in the landscape
  nodes_remaining_by_cell = list() # a list of cells with all the genes remaining in each cell
  time=0 # backward time
  single_coalescence_events=0 # number of single coalescence events. Coalescence involving multiple individuals counts for 1 event.
  single_and_multiple_coalescence_events=0 # number of single and multiple coalescence events. Coalescence involving multiple individuals counts for "the number of individuals - 1" events.
  for (cell in 1:ncell(rasterStack))#cell=1)
  {
    nodes_remaining_by_cell[[cell]] <- which(cell_number_of_nodes==cell)
  }
  while (length(unlist(nodes_remaining_by_cell))>1) 
  {
    # migration
    # we localize the parents in the landscape by sampling in the backward transition matrix
    for (node in 1:length(parent_cell_number_of_nodes))#gene=1;node=1# parent_cell_number_of_nodes
    {
      parent_cell_number_of_nodes[node] = sample(ncell(rasterStack),size=1,prob=c(transitionmatrice[cell_number_of_nodes[node],]))
    }
    # once we know the parent cell numbers, we calculate the forward dispersion probability of the event
    prob_forward[time] = sum(log( transition_forward[parent_cell_number_of_nodes,cell_number_of_nodes]))
    # coalescence
    time=time+1; if (round(time/10)*10==time) {print(time)}
    # we now perform coalescence within each cell of the landscape for the parents
    for (cell in 1:ncell(rasterStack))#cell=1;cell=2;cell=3;cell=4;cell=5;cell=26;cell=10
    {     
      nodes_remaining_in_the_cell = nodes_remaining_by_cell[[cell]] <- as.numeric(names(which(parent_cell_number_of_nodes==cell)))
      # we obtain the identities in the geneticData table (line) of the genes remaining in the cell
      if (length(nodes_remaining_in_the_cell)>1) 
      {
        # We create a logical matrix in which lines represent genes of the sample (nodes) remaining in the cell 
        # and column represent their parent chosen from the whole population of size N[cell]. 
        # If two genes (lines) coalesce if they have TRUE for the same parent (column) 
        nbgenesremaining=length(nodes_remaining_in_the_cell)
        smp = sample(N[cell],length(nodes_remaining_in_the_cell),replace=TRUE)
        parentoffspringmatrix <- matrix(smp,nrow=nbgenesremaining,ncol=N[cell])==matrix(1:N[cell],nrow=nbgenesremaining,ncol=N[cell],byrow=TRUE)
        #        colnames(parentoffspringmatrix) <- nodes_remaining_in_the_cell
        rownames(parentoffspringmatrix) <- nodes_remaining_in_the_cell
        # the columns  column of parentoffspringmatrix that have more than one TRUE
        # identifies the individuals that coalesce with the lines of the TRUEs
        if (any(colSums(parentoffspringmatrix)>1) )
        {
          #  which(colSums(parentoffspringmatrix)>1)) gives the column names 
          #  of parentoffspringmatrix that have coalescence information
          for (multiple in which(colSums(parentoffspringmatrix)>1)) # multiple<-which(colSums(parentoffspringmatrix)>1)[1]
          {
            # there is coalescence 
            single_coalescence_events = single_coalescence_events +1
            # which(parentoffspringmatrix[,multiple]) identifies which node in the column coalesce
            nodes_that_coalesce = names(which(parentoffspringmatrix[,multiple]))
            # attibutes new node number to the ancestor, adds this to the nodes vector, removes the nodes that coalesced from the node vector
            new_node <- max(nodes)+1;nodes=nodes[!(names(nodes)%in%nodes_that_coalesce)];nodes=append(nodes,new_node);names(nodes)[length(nodes)]=new_node
            # updating of vector parent_cell_number_of_nodes (adding the cell number of the new node and removing the nodes that disapeared)
            parent_cell_number_of_nodes <- append(parent_cell_number_of_nodes[!(names(parent_cell_number_of_nodes)%in%nodes_that_coalesce)],cell);names(parent_cell_number_of_nodes)[length(parent_cell_number_of_nodes)]<-new_node
            # adds the event to the list coalescent: time, which node coalesced, and the number of the new node
            coalescent[[single_coalescence_events]] <- list(time=time,coalescing=as.numeric(nodes_that_coalesce),new_node=new_node)
            # updating the nodes vector for the cell
            nodes_remaining_in_the_cell = nodes_remaining_by_cell[[cell]] <- append(nodes_remaining_in_the_cell[!nodes_remaining_in_the_cell %in% nodes_that_coalesce],new_node)
            # updates the number of coalescent events 
            single_and_multiple_coalescence_events = single_and_multiple_coalescence_events + length(nodes_that_coalesce) - 1
          }
        }
      }
    }
  # we now move in the backward generation while coalescence loop
  cell_number_of_nodes = parent_cell_number_of_nodes
  }
  list(coalescent=add_br_length_and_mutation(coalescent,mutation_rate),forward_log_prob=sum(prob_forward)/coalescent[[length(coalescent)]]$time)
  # forward_log_prob is the average per generation of the log probability of the forward movements of the genes in the landscape
}

#
# add br_length and mutation to coalescent list
#
#

add_br_length_and_mutation <- function(coalescent,mutation_rate=1E-4)
{
  tips = NULL
  internals = NULL
  nodes = NULL
  times = NULL
  for (i in 1:length(coalescent))#i=1;i=2
  {
    nodes = append(nodes,coalescent[[i]]$coalescing,coalescent[[i]]$new_node)
    internals = append(internals,coalescent[[i]]$new_node)
    times = append(times,coalescent[[i]]$time)
  }
  nodes = as.numeric(levels(as.factor(c(nodes,internals))));nodes = nodes[order(nodes)]
  tips = nodes[!((nodes)%in%(internals))]
  # getting the branch length of each coalescing node
  for (i in 1:length(coalescent))#i=1
  {
    for (coalescing in coalescent[[i]]$coalescing)# coalescing = coalescent[[i]]$coalescing[1]
    {
      if (coalescing %in% tips) {coalescent[[i]]$br_length <- append(coalescent[[i]]$br_length,coalescent[[i]]$time)
                                                                    } else {
        coalescent[[i]]$br_length <- append(coalescent[[i]]$br_length,coalescent[[i]]$time-times[which(internals==coalescing)]) 
                                                                    } 
      coalescent[[i]]$mutations <- rpois(rep(1,length(coalescent[[i]]$br_length)),coalescent[[i]]$br_length*mutation_rate)
    }
  }  
coalescent
}

# coalescent_2_newick
# fuinction that converts coalescent
# 

coalescent_2_newick <- function(coalescent)
{
  #  tree=paste("(",coalescent[[length(coalescent)]]$new_node,")",sep="")
  tree=paste(" ",coalescent[[length(coalescent)]]$new_node," ",sep="")
  for (i in length(coalescent):1)
  {
    Time = coalescent[[i]]$time
    coalesc <- as.character(coalescent[[i]]$coalescing)
    tree <- str_replace(tree,paste(" ",as.character(coalescent[[i]]$new_node)," ",sep=""),paste(" ( ",paste(" ",coalesc," :",coalescent[[i]]$br_length,collapse=" ,",sep=""),") ",sep=""))
  }
  tree <- gsub(" ","",paste(tree,";",sep=""))
tree
}

# coalescent_2_phylog
# converts a coalescent simulation to a class phylog tree (library ape)
#

coalescent_2_phylog <- function(coalescent)
{
read.tree(text=coalescent_2_newick(coalescent))
}

#
# plot_coalescent plots a coalescent simulation
# argument: output of simul_coalescent()

plot_coalescent <- function(coalescent,with_landscape=FALSE,rasK=NULL)
{
  if (with_landscape) {par(mfrow=c(1,2),oma=c(0,0,0,4),xpd=TRUE)}else{par(mfrow=c(1,1),oma=c(0,0,0,4),xpd=TRUE)}
  tipcells <- geneticData$Cell_numbers[as.numeric(coalescent_2_phylog(coalescent)$tip.label)]
  tipcols = rainbow(ncell(rasK))[tipcells]
  plot(coalescent_2_phylog(coalescent),direction="downward",tip.color=tipcols)
  legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(-1,0))
  if (with_landscape) {plot(rasK)}
}

#
# genetic_simulation
# this function simulates genetic data from a coalescent
#
#


test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice){
  ## ref: Rousset et al. J Evol Biol,13  (2000) 58-62.
  vecteur_a_value <-c(0)
  GenDist=list()
  for(i in 1: nb_generations){#i=1
    print(i)
    geneticData <- repnDispMutFunction(geneticData, dimGeneticData, mutationRate, transitionmatrice) 
    Genotypes = geneticData[,grep("Locus", colnames(geneticData), fixed = T)]
    popmbrship=geneticData[,"Cell_numbers"]
    Fst = Fstat(Genotypes,nCell,popmbrship,ploidy=2)
    Fis = Fst[[1]];Fst=Fst[[2]]
    GenDist[[i]] <- Fst/(1-Fst)
    if((i>90) && (i%%30 == 0)){
      if(var(GenDist[[(i-30):i]])> var(GenDist[[(i-60):(i-30)]]) 
         && var(GenDist[[(i-30):i]])> var(GenDist[[(i-90):(i-60)]])){
        return(list(geneticData, GenDist))
               break
      }
    }
  }
}

# to get a matrix of a-values between all cells of a landscape from 
# simulations of evolution where not all cells necesarily contain individuals
#
#

a_value_matrix_from_forward_simulation <- function(geneticFinal,rasterStack)
{
  matrixQb = (1-Qwithin_pair(geneticFinal)+2*(Qwithin_pair(geneticFinal)-Qbetween(geneticFinal,dim(geneticFinal))))
  matrixQw = 2*(1-Qwithin_pop(geneticFinal))
  vecteur_a_value = matrixQb/matrixQw-1/2
  vecteur_a_value[is.na(vecteur_a_value)] <-0
  cell_numbers_list <- levels(as.factor(geneticFinal$Cell_numbers))
  agg <- aggregate(vecteur_a_value,by=list(geneticFinal$Cell_numbers),FUN="mean")[,-1]
  agg2 <- t(aggregate(t(agg),by=list(geneticFinal$Cell_numbers),FUN="mean"))[-1,]
  row.names(agg2) <- as.numeric(cell_numbers_list)
  colnames(agg2) <- as.numeric(cell_numbers_list)
  agg3 <- matrix(NA, nrow=ncell(rasterStack),ncol=ncell(rasterStack))
  agg3[as.numeric(row.names(agg2)),as.numeric(colnames(agg2))] <- agg2
agg3
}

# Function that determine number of dispersal parameters from dispersal shapeDisp used
# Useful for nlm estimation
nbpDisp <- function(shapeDisp){
  (switch(shapeDisp,
          fat_tail1 = 2,
          gaussian = 1,
          exponential = 1,
          contiguous = 1,
          island = 1,
          fat_tail2 = 2))
}

# Function that simulate a genetic data with parameters given. 
#It returns a list of 2 variables:  final genetic data observed and matrix of a-value observed
simulationGenet <- function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, mutationRate, nbLocus, initial_locus_value, shapeDisp, pDisp, nb_generations,ind_per_cell=30){
  K <- subset(ReactNorm(values(donneesEnvironmentObs), pK , shapesK),select="Y") # carrying capacity
  r <- subset(ReactNorm(values(donneesEnvironmentObs), pR , shapesR),select="Y") # growth rate
  Rast_K <- donneesEnvironmentObs ; values(Rast_K) <- K
  Rast_r <- donneesEnvironmentObs ; values(Rast_r) <- r
  geneticData = CreateGenetArray(Rast_K, nbLocus,initial_locus_value,Option="full_pop")
  geneticData[,"Cell_number_init"] <- geneticData[,"Cell_numbers"]
  dimGeneticData = dim(geneticData)
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
  transitionmatrice = transitionMatrixBackward(Npop = K, migration= migrationM)
  geneticDataFinal = test_stabilite_a_value(geneticData, mutationRate, dimGeneticData, nb_generations,transitionmatrice)
  a_value_obs = geneticDataFinal[[2]]
  geneticDataFinal = geneticDataFinal[[1]]
  return(list(geneticDataFinal, a_value_obs))
}

# plots matrixes of forward mutation for validation
#
#

matrixes_forward <- function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, shapeDisp, pDisp, a_value_obs, a_value_att, file=NULL, mutationRate,nbLocus, initial_locus_value,indpercell)
{
  nblayers =dim(donneesEnvironmentObs)[3]
  nCell = ncell(donneesEnvironmentObs)
  Cell_numbers <- 1:nCell
  geneticObs = simulationGenet(donneesEnvironmentObs,pK,pR,shapesK,shapesR,mutationRate,nbLocus,initial_locus_value,shapeDisp,pDisp,nb_generations=5000,indpercel)
  finalGenetData = geneticObs[[1]]
  a_value_simul = geneticObs[[2]]
  K <-reactNorm(donneesEnvironmentObs, pK, shapesK) # carrying capacity
  r <-reactNorm(donneesEnvironmentObs, pR, shapesR) # carrying capacity
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
  transitionmatrice = transitionMatrixForward(r, K, migration= migrationM)
  land_size <- raster(matrix(K[1,],nrow=dim(donneesEnvironmentObs)[1],ncol=dim(donneesEnvironmentObs)[2]))# be careful transposed
  land_size <- raster(matrix(K[1,]))
  extent(land_size) <- extent(donneesEnvironmentObs)
  list(land_size,transitionmatrice,a_value_simul,a_value_att)
}

# MAIN DES FONCTIONS POUR SIMULER DES DONNEES GENETIQUES OBSERVEES
#mutationRate = 10^(-4);                             shapeDisp = "fat_tail1";
#nbLocus=10;                                         initial_locus_value=200
#pDisp = c(alphaDisp = 20, betaDisp = 1.5);     nbpDisp =nbpDisp(shapeDisp);
#aggregate_factor=16
#donneesEnvironmentObs = Aggregate_and_adjust_raster_to_data(raster(paste(wd,envdir,envfiles,sep="")),release=read.table(paste(wd,genetfile,sep="")), recovery=read.table(paste(wd,genetfile,sep="")), extend_band_size=1, aggregate_index=aggregate_factor)
  
# donneesEnvironmentObs = raster(matrix(NA,nrow=20,ncol=20))
#values(donneesEnvironmentObs) <-log(ceiling(runif(ncell(donneesEnvironmentObs),0,1)*3))
#plot(donneesEnvironmentObs)
#nblayers =dim(donneesEnvironmentObs)[3]
#nCell = ncell(donneesEnvironmentObs)
#geneticObs = simulationGenet(donneesEnvironmentObs,alpha=0, beta =1,mutationRate,nbLocus, initial_locus_value,shapeDisp,pDisp,nb_generations=5000)
#finalGenetData = geneticObs[[1]]
#a_value_obs = geneticObs[[2]]

#####################################
##### ESTIMATION DES PARAMETRES #####
#####################################
expect_a_value <- function(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,Cell_numbers,nbpDisp,nblayers,shapeDisp)
{
  K = ReactNorm(rasterStack =donneesEnvironmentObs, pK, shapesK)
  r = ReactNorm(rasterStack =donneesEnvironmentObs, pR, shapesR)
  matrice_migration = migrationMatrix(donneesEnvironmentObs,shapeDisp, shapeDisp)
  matrice_transition = transitionMatrixBackward(K,r, matrice_migration)
  matrice_laplace = laplaceMatrix(matrice_transition)
  commute_time = resistDist(matrice_laplace)
  genetic_dist_att = geneticDist(commute_time, popSize)
  a_value_att = genetic_dist_att[Cell_numbers,Cell_numbers]
  a_value_att
}

ssr <- function(p){
  a_value_att = expect_a_value(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,Cell_numbers=finalGenetData$Cell_Number,nbpDisp,nblayers,shapeDisp)
  return(mean((a_value_obs - a_value_att)^2))
}

#initial = c(0,1,20,1.5)
#fct_erreur_calc = nlm(f = ssr , p = initial, hessian = FALSE)
#p=initial
#a_value_graph_model = expect_a_value(donneesEnvironmentObs,initial,finalGenetData$Cell_Number,nbpDisp=2,nblayers=1,shapeDisp="fat_tail1")
#ssr(initial)
#parametres_infere = fct_erreur_calc$estimate
#parametres_reels = c(alpha = 0, beta = 2,beta = 1, alphaDisp = 20, betaDisp = 1.5)

validation <- function(donneesEnvironmentObs,pK = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))), pR = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))),shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"),shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"), shapeDisp="fat_tail1", pDisp = c(mean=0.32,shape=1.6), file=NULL,mutationRate,nbLocus, initial_locus_value,nb_generations=5000,indpercell=30)
{
  nblayers =dim(donneesEnvironmentObs)[3]
  nCell = ncell(donneesEnvironmentObs)
  Cell_numbers <- 1:nCell
  K <- ReactNorm(values(donneesEnvironmentObs), pK , shapesK)$Y # recuperation de la carrying capacity
  r <- ReactNorm(values(donneesEnvironmentObs), pR , shapesR)$Y # recuperation du growth rate
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
  transitionmatrice = transitionMatrixBackward(Npop = K, migration= migrationM)
  geneticObs = simulationGenet(donneesEnvironmentObs,alpha, beta, mutationRate,nbLocus, initial_locus_value,shapeDisp,pDisp,nb_generations=5000,indpercell)
  finalGenetData = geneticObs[[1]]
  land_size <- raster(matrix(K[1,],nrow=dim(donneesEnvironmentObs)[1],ncol=dim(donneesEnvironmentObs)[2]))
  a_value_simul = a_value_matrix_from_forward_simulation(geneticObs[[1]],land_size)
  a_value_theory_stepping_stone_model = distanceMatrix(donneesEnvironmentObs)/(4*0.05)
  a_value_theory_island_model <- matrix(pDisp[1],nrow=nCell,ncol=nCell)-pDisp[1]*diag(nCell)
  a_value_theory_graph_model <- expect_a_value(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,Cell_numbers,nbpDisp=nbpDisp,nblayers=nblayers,shapeDisp)
  if (!is.null(file)){jpeg(file)}
  par(mfrow=c(2,3))
  plot(land_size,main="Population sizes")
  #plot(raster(migrationM))
  plot(raster(transitionmatrice),main="Transition matrix")
  plot(raster(a_value_simul),main="Simul genet differ")
  plot(raster(a_value_theory_stepping_stone_model),main="Expected stepping stone")
  plot(raster(a_value_theory_island_model),main="Expected island")
  plot(raster(a_value_theory_graph_model),main="Expected graph model")
  if (!is.null(file)){dev.off()}
  list(land_size,transitionmatrice,a_value_simul,a_value_theory_stepping_stone_model,a_value_theory_island_model,a_value_theory_graph_model,finalGenetData)
}


#test nlm
ssr <- function(p){
  popSize = K_Function(rasterStack = donneesEnvironmentObs, p[1], p[2:(nblayers+1)])
  matrice_migration = migrationMatrix(donneesEnvironmentObs,shapeDisp, p[(nblayers+2):(nbpDisp+nblayers+1)])
  matrice_transition = transitionMatrixBackward(popSize, matrice_migration)
  matrice_laplace = laplaceMatrix(matrice_transition)
  commute_time = resistDist(matrice_laplace)
  genetic_dist_att = geneticDist(commute_time, p[(nbpDisp+nblayers+2)])
  
  a_value_att = genetic_dist_att[finalGenetData$Cell_Number,finalGenetData$Cell_Number]
  return(mean((a_value_obs - a_value_att)^2))
}



