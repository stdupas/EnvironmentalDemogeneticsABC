#
##########################################################################

# matrice : colonne
#raster : ligne


#############################
#### Library to download ####
#############################

library(rgdal)
library(raster)
library(MASS)
#library(Geneland)
library(ape)
library(stringr)
library(lattice)
library(markovchain)
library(matrixcalc)
library(abind)


####################################################
##### BACKWARD MODEL FUNCTIONS AND EXECUTIONS ######
#####  SIMULATION OF PREDICTED GENETIC DATA   ######
####################################################


# Function to get spatial resolution in km
degree2km <- function(rasterStack)
{
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

proportional <- function(X,p,Log=FALSE)
{
  if (Log) {log(p[rep("a",dim(X)[1]),colnames(X)]*X)} else {
    p[rep("a",dim(X)[1]),colnames(X)]*X
  }
}








################################################################################################
#########################     METHODE A IMPLEMENTER     ########################################
################################################################################################














linear <- function(X,p,Log=FALSE)
{
  Yx1 = p[rep("Yx1",dim(X)[1]),colnames(X)]
  Yx0 = p[rep("Yx0",dim(X)[1]),colnames(X)]
  a = (Yx1 - Yx0)
  b = Yx0
  if (Log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
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
  # X is a matrix with the different environemental variables in column and the map cells in line
  # or a numeric vector
  #p=c(p["Xmin"]=10,p["Xmax"]=20,p["Xopt"]=18,p["Ymax"]=0.1)
  ## shapeDisps in c("enveloppe","envelin","envloglin","loG","conquadratic","conquadraticskewed","conquadraticsq","conquadraticskewedsq")
{
  if (class(p)=="numeric") {p=t(as.matrix(p));rownames(p)="a"}
  if (class(X)=="numeric") {X=data.frame(X=X);colnames(X)=colnames(p)}
  Y=X
  if (!all(colnames(p)%in%names(shapes))) {stop ("variable names do not correspond between parameters 'p' ans 'shape'")}
  for (shape in as.character(levels(as.factor(shapes))))
  {
    variables = colnames(p)[which(shapes==shape)]
    Y[,variables]=switch(shape,
                         constant=p,
                         proportional = proportional(subset(X,select=variables),p),
                         linear = linear(subset(X,select=variables),p),
                         enveloppe=enveloppe(subset(X,select=variables),p),
                         envelin=envelinear(subset(X,select=variables),p),
                         envloglin=envelinear(subset(X,select=variables),p,log=TRUE),
                         loG = log(subset(X,select=variables)),
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
K_Function <- function(rasterStack, p, shapes)
{
  #K=matrix(combineReactNorms(values(rasterStack),p,shapes),byrow=TRUE,nrow=length(rasterStack),ncol=legnth(rasterStack))
  ReactNorm(valuesA(rasterStack),p,shapes)
}

# R_Function: Gets effective growth rate from environmental variables 
# alpha and beta are fixed by estimation
# rasterStack = environmental variables 
R_Function <- function(rasterStack, alpha, beta)
{
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
  populationSize <- donneesEnvironmentObs[[1]]
  populationSize[cellNumA(populationSize)] <- ReactNorm(valuesA(donneesEnvironmentObs), p, shapes)[1,]
  populationSize
}
  
# (optional) distanceMatrix return distance between all cells of raster
distanceMatrix <- function(rasterStack)
{
  #get x and y coordinates for each cell of raster object put in parameters
  coords = xyFromCellA(rasterStack)
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
migrationMatrix <- function(rasterStack,shapeDisp, pDisp)
{
  distMat<-distanceMatrix(rasterStack)
  Ndim = 1+all(ncell(rasterStack)!=dim(rasterStack)[1:2])
  migration = apply(distMat, c(1,2), 
                    function(x)(switch(shapeDisp,
                                       # 1: alphaDisp   2: betaDisp ; note: esperance = 1/alphaDisp
                                       fat_tail1 = 1/(1+x^pDisp[2]/pDisp[1]), # Molainen et al 2004
                                       # 1: sigmaDisp
                                       gaussian = (dnorm(x, mean = 0, sd = pDisp[1], log = FALSE)),
                                       # 1: sigma Disp
                                       exponential = (dexp(x, rate = 1/pDisp[1], log = FALSE)),
                                       # 1: sigmaDisp
                                       contiguous = (x==0)*(1-pDisp[1])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp[1]/(2*Ndim)),
                                       # 1: sigmaDisp
                                       contiguous8 = (x==0)*(1-pDisp[1])+((x>0)-(x>2*res(rasterStack)[1]))*(pDisp[1]/(4*Ndim)),
                                       #island = (x==0)*(1-pDisp[1])+(x>0)*(pDisp[1]/(nCell-1)),
                                       island = (x==0)*(1-pDisp[1])+(x>0)*(pDisp[1]),
                                       #1: sigmaDisp    2: gammaDisp
                                       fat_tail2 = x^pDisp[2]*exp(-2*x/(pDisp[1]^0.5)),
                                       #1: sigmaDisp, 2: 
                                       contiguous_long_dist_mixt = pDisp["plongdist"]/ncellA(rasterStack)+(x==0)*(1-pDisp["pcontiguous"]-pDisp["plongdist"])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp["pcontiguous"]/2),
                                       gaussian_long_dist_mixt = pDisp[2]/ncellA(rasterStack) + (dnorm(x, mean = 0, sd = pDisp[1], log = FALSE))
                    )))
  return(migration)
}

# transitionMatrix obtained with an isotropic migration hypothesis for a backward model
transitionMatrixBackward <- function(r,K, migration)
{
  if ((length(r)==1)&(length(K)==1)){transition = r * K * t(migration)}
  if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
  if ((length(r)==1)&(length(K)>1)){transition = r * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
  if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
  if ((length(r)>1)&(length(K)>1)) {transition = t(matrix(r,nrow=length(r),ncol=length(r))) * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
  transition = transition / rowSums(transition)  # Standardisation
  transition = transition * as.numeric(transition>1e-6) # removal of values below 1e-6
  transition = transition / rowSums(transition)  # Standardisation again
  transition
}

# transitionMatrix obtained with an isotropic migration hypothesis for a formard model
transitionMatrixForward <- function(r,K, migration, meth="non_overlap")
{
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

#
# Absorbing transition matrix
# Estimates probability of transition among transient and absorbant states
# for pairs of genes in demes newtork
# dimension is N(N+3)/2
# arg: 
# transition : backward transition probability matrix of genes among demes
# N : population sizes of the demes
# Mcrae 2006
# Hey 1001

absorbingTransition <- function(transition,N)
{
  N[N<1]=1
  Ndeme <- dim(transition)[1]
  # States of pairs of genes
  # N Homodemic not coalesced states
  # Heterodemic states (ij)!=(kl)
  Nhetero <- Ndeme*(Ndeme-1)/2
  # Ndeme Coalesced demic state, can be merged to 1 if we don't mind the deme
  #
  # First: 
  # Calculate the transient states (not coalesced) transition matrix
  # for transition among Q among N*(N+1)/2 not coalesced states 
  # It is composed of a submatrixes of transition between not coalesced heterodemic and
  # homodemic states
  #
  #    /          \
  #    |HeHe  HeHo|
  # Q= |HoHe  HoHo|
  #    \          /
  # where He are heterodemic states, and Ho are homodemic states
  # In submatrix HeHe:
  # lines ij are ordered for (i in 2:Ndeme) for (j in 1:(i-1))
  # columns kl are ordered for (k in 2:Ndeme) for (l in 1:(k-1))
  # In submatrix HoHe the lines are from 1 to Ndeme
  Qline <- matrix(NA,Ndeme,Ndeme) # this matrix of the same dimension as the transition matrix
                                  # contains the line number or column number 
                                  # in transient Q matrix for each pair of deme {i,j} 
                                  # presented as in the transition matrix
                                  # 
  QheteroHetero <- matrix(NA,Nhetero,Nhetero)
  ij=0;kl=0
  Check=TRUE
  for (i in 2:Ndeme){
    for(j in 1:(i-1)) {
      ij=ij+1
      Qline[i,j] <- Qline[j,i] <- ij # this matrix aims to find where the pair of demes {i,j}
                                     # is in Q matrix lines or columns
#      i_j_Q_lines[ij,]<-c(i,j)
    kl=0
    for (k in 2:Ndeme){
        for (l in 1:(k-1)){
          kl=kl+1
          QheteroHetero[ij,kl] <- transition[i,k]*transition[j,l]
#          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
        }
      }
    }
  }

  # then transient matrix 
  QhomoHetero <- matrix(NA,Ndeme,Nhetero)
  kl=0
  Check=TRUE
  for (i in 1:Ndeme){ # only homodemic sources are considered (i=j)
    Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
                            # in the lines of Q matrix
    for (k in 2:Ndeme){
      for (l in 1:(k-1)){ # only heterodemic targets
        kl=kl+1
        QhomoHetero[i,kl] <- transition[i,k]*transition[i,l]    # i=j (homodemic sources)
        #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
      }
    }
    kl=0
  }
  QheteroHomo <- matrix(NA,Nhetero,Ndeme)
    ij=0
    for (i in 2:Ndeme){
      for (j in 1:(i-1)){ # only heterodemic sources
        ij=ij+1
        for(k in 1:Ndeme){ # only homodemic targets
          QheteroHomo[ij,k] <- transition[i,k]*transition[j,k]*(1-1/(2*N[k]))
          # homodemic targets that have not coalesced
      }
    }
  }
  QhomoHomo <- transition*transition*matrix(1-1/(2*N),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
Q <- cbind(rbind(QheteroHetero,QhomoHetero),rbind(QheteroHomo,QhomoHomo))
list(Q=Q,Qline=Qline)
}

#
#
#
#

absorbingTransitionComplete <- function(transition,N)
{
  N[N<1]=1
  Ndeme <- dim(transition)[1]
  # States of pairs of genes
  # N Homodemic not coalesced states
  # Heterodemic states (ij)!=(kl)
  Nhetero <- Ndeme*(Ndeme-1)/2
  # Ndeme Coalesced demic state, can be merged to 1 if we don't mind the deme
  #
  # First: 
  # Calculate the full transient states (not coalesced) and absorbant states (coaleced) 
  # transition matrix
  # for transition among Q among N*(N+1)/2 not coalesced states 
  # It is composed of a submatrixes of transition between not coalesced heterodemic and
  # homodemic states
  # 
  #    /              \
  #    |HeHe HeHo HeCo|
  #    |HoHe HoHo HoCo|
  #    |CoHe CoHo Coco|
  #    \              /
  # 
  #
  #
  # where He are heterodemic states, and Ho are homodemic states
  # In submatrix HeHe:
  # lines ij are ordered for (i in 2:Ndeme) for (j in 1:(i-1))
  # columns kl are ordered for (k in 2:Ndeme) for (l in 1:(k-1))
  # In submatrix HoHe the lines are from 1 to Ndeme
  Qline <- matrix(NA,Ndeme,Ndeme) # this matrix of the same dimension as the transition matrix
  # contains the line number or column number 
  # in transient Q matrix for each pair of deme {i,j} 
  # presented as in the transition matrix
  # 
  HoHe <- matrix(NA,Ndeme,Nhetero)
  HeHo <- matrix(NA,Nhetero,Ndeme)
  CoHe <- matrix(0,Ndeme,Nhetero)
  HeCo <- matrix(NA,Ndeme,Nhetero)
  HoCo <- matrix(transition/(2*matrix(N,nrow=Ndeme,ncol=Ndeme, byrow=TRUE)),nrow=Ndeme,ncol=Ndeme)
  CoHo <- matrix(0,nrow=Ndeme,ncol=Ndeme)
  HeHe <- matrix(NA,Nhetero,Nhetero)
  CoCo <- transition
  HoHo <- transition*transition
  ij=0;kl=0
  Check=TRUE
  for (i in 2:Ndeme){
    for(j in 1:(i-1)) {
      ij=ij+1
      Qline[i,j] <- Qline[j,i] <- ij # this matrix aims to find where the pair of demes {i,j}
      # is in Q matrix lines or columns
      #      i_j_Q_lines[ij,]<-c(i,j)
      kl=0
      for (k in 2:Ndeme){
        for (l in 1:(k-1)){
          kl=kl+1
          HeHe[ij,kl] <- transition[i,k]*transition[j,l]
          #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
        }
      }
      for (k in 1:Ndeme){
        HeHo[ij,k] <- transition[i,k]*transition[j,k]*(1-1/(2*N[k]))
        HeCo[ij,k] <- transition[i,k]*transition[j,k]/(2*N[k])
        HoHe[k,ij] <- transition[k,i]*transition[k,j]
      }
     }
    }
  
  # then transient matrix 
  kl=0
  Check=TRUE
  for (i in 1:Ndeme){ # only homodemic sources are considered (i=j)
    Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
    # in the lines of Q matrix
    for (k in 2:Ndeme){
      for (l in 1:(k-1)){ # only heterodemic targets
        kl=kl+1
        HoHe[i,kl] <- transition[i,k]*transition[i,l]    # i=j (homodemic sources)
        #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
      }
    }
    kl=0
  }
  QheteroHomo <- matrix(NA,Nhetero,Ndeme)
  HeCo <- matrix(NA,nrow=Nhetero,ncol=Ndeme)
  ij=0
  for (i in 2:Ndeme){
    for (j in 1:(i-1)){ # only heterodemic sources
      ij=ij+1
      for(k in 1:Ndeme){ # only homodemic targets
        QheteroHomo[ij,k] <- transition[i,k]*transition[j,k]*(1-1/(2*N[k]))
        HeCo[ij,k] <- transition[i,k]*transition[j,k]/(2*N[k])
        # homodemic targets that have not coalesced
      }
    }
  }
  QhomoHomo <- transition*transition*matrix(1-1/(2*N),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
  Q <- cbind(rbind(QheteroHetero,QhomoHetero),rbind(QheteroHomo,QhomoHomo))
  CoHeCoHo <- matrix(0,nrow=Ndeme,ncol=Nhetero+Ndeme)
  #HeCoHoCo <- matrix(NA,nrow=Nhetero+Ndeme,ncol=Ndeme)
  ij=0
#  for (i in 2:Ndeme){
#    for(j in 1:(i-1)) {
#      ij=ij+1
#      for (k in 1:Ndeme)
#      }
#  }
  list(Q=Q,Qline=Qline)
}

#
# Fundamental matrix of absorbing chain
#

fundamentalMatrixAbsorbingChain <- function(transientTransitionMatrix)
{# transtientTransitionMatrix
  Ndeme <- (2*dim(transientTransitionMatrix)[1]+1/4)^.5-.5
  solve(diag(Ndeme*(Ndeme+1)/2)-transientTransitionMatrix)
}

#
# Absorbing times vector and matrix
#

absorbingTime <- function(fundamentalMatrixAbsChain,Qline)
{
  # first use the fundamental matrix to calculate
  # absorbingtimes for pair of genes {i,j} in the following sequence:
  # 1) for (i in 2:Ndeme) for (j in 1:(i-1)) (Heterodemic states)
  # 2) for (i=j in 1: Ndeme) (homodemic states)
  Ndeme <- (2*dim(fundamentalMatrixAbsChain)[1]+1/4)^.5-.5
  absTimesVector <- fundamentalMatrixAbsChain%*%rep(1,nrow(fundamentalMatrixAbsChain))
  # then put these absorbing times in a matrix Ndeme by Ndeme
  absTimesMatrix <- matrix(absTimesVector[Qline],nrow=Ndeme)
  absTimesMatrix
}

#
# genetic Distance
#

genetDistAbsorbingMethod <- function(transition,N,mutation_rate)
{
  # N is effective population size of the deme
  # transition is the backward transition matrix among demes
  Qlist<-absorbingTransition(transition,N)
  absorbingTime(fundamentalMatrixAbsorbingChain(Qlist$Q),Qlist$Qline)*mutation_rate
}

#
#
#

coalescenceTimeProbDistrib <- function(Qlist)
{
  x <- eigen(Qlist[["Q"]])
  Lambda <- x$values;Theta=x$vector
  Ndeme <- (2*dim(Q)[1]+1/4)^.5-.5
  cTPD <- array(NA,dim=c(Ndeme,Ndeme,1))
  condition=TRUE
  t=0
  LambdaPower <- rep(1,length(x$values))
  while (condition){
    t=t+1
    LambdaPowerBefore <- LambdaPower
    LambdaPower <- LambdaPowerBefore*Lambda
    cTPD[,,t] <- matrix(colSums(Re(Theta%*%diag(LambdaPowerBefore-LambdaPower)%*%t(Theta))[,c(Qlist[["Qline"]])]),Ndeme,Ndeme)
  }
}

#
# Probability distribuiton of coalescence times
# Hey 2001

coalescenceTimeProbDistrib <- function(Qlist)
{
  Ndeme <- (2*dim(Qlist[[1]])[1]+1/4)^.5-.5
  cTPD <- array(NA,dim=c(Ndeme,Ndeme,1))
  condition=TRUE
  t=1
  QPower <- Q <- Qlist[["Q"]]
  cTPD[,,t] <- matrix(colSums((matrix(1,nrow(QPower),ncol(QPower))-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme)
  while (condition){
    t=t+1
    QPowerBefore <- QPower
    QPower <- Q%*%QPowerBefore
    cTPD <- abind(cTPD,matrix(colSums((QPowerBefore-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme),along=3)
    condition=!all(colSums(cTPD,3)>0.9999)
  }
cTPD
}

#
# Probability of coalescence time
# Hey 2001

coalescenceProb <- function(Qlist,time)
{
  Ndeme <- (2*dim(Qlist[[1]])[1]+1/4)^.5-.5
  cTPD <- array(NA,dim=c(Ndeme,Ndeme,1))
  condition=TRUE
  t=1
  QPower <- Q <- Qlist[["Q"]]
  cTPD[,,t] <- matrix(colSums((matrix(1,nrow(QPower),ncol(QPower))-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme)
  while (condition){
    t=t+1
    QPowerBefore <- QPower
    QPower <- Q%*%QPowerBefore
    cTPD <- abind(cTPD,matrix(colSums((QPowerBefore-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme),along=3)
    condition=!all(colSums(cTPD,3)>0.9999)
  }
  cTPD
}

# simulation forward of population sizes across time
forward_simul_landpopsize <- function(N0,p, migration)
{
  
}
# laplaceMatrix returns Laplacian matrix from transition matrix

# DONE !!!!!!!!!!!!!!!!!!
laplaceMatrix <- function(transitionMatrix)
{
  matrixD = diag(rep(1,dim(transitionMatrix)[1])) # diagonal equals to 1
  laplacianMatrix = matrixD - transitionMatrix
  laplacianMatrix[is.na(laplacianMatrix)]<-0 # replace NA by 0
  laplacianMatrix
}
# ordinary laplacian
# Boley et al 2011
#

# DONE !!!!!!!!!!!!!!!!
ordinary_laplacian <- function(transition)
{
  markovB<-new("markovchain", states=dimnames(transition)[[1]], transitionMatrix=transition)
  PI<-diag(steadyStates(markovB)[1,])
  PI - PI%*%transition
}

# Calcul of resistance between two points of the graph
# with the Moore-Penrose generalized inverser matrix.
# ref : Bapat et all, A Simple Method for Computing Resistance Distance (2003)
# ref : Courrieu, Fast Computation of Moore-Penrose Inverse Matrices (2005)
#r <- inverseMP
#for (i in 1:dim(r)[1]){
#  for (j in 1:dim(r)[2]){
#    r[i,j] <- inverseMP[i,i]+inverseMP[j,j]-inverseMP[i,j]-inverseMP[j,i]
#  }
#}
  

# DONE !!!!!!!!!!!!!!
commute_time_undigraph <- function(matrice_transition)
{
  laplacian = laplaceMatrix(matrice_transition)
  inverseMP = ginv(laplacian) # generalized inverse matrix  (Moore Penrose)
  diag = diag(inverseMP) # get diagonal of the inverse matrix
  mii = matrix(diag, nrow =dim(inverseMP), ncol = dim(inverseMP))
  mjj = t(mii)
  mij = inverseMP
  mji = t(mij)
  commute_time = mii + mjj - mij - mji
  commute_time
}

# for weighted digraphs
# Boley et al 2011
# Hitting time is given by 
# H = 1 · [diag(Z )]T − Z 
# Where Z, the scaled Fundamental Matrix is defined by
# Z = (L + pi_ pi_ T)⁻¹ = (PI ( 1 - P π T))⁻¹
# where P is the transition probability matrix, 
# pi_ (or π) is the stationary probablity distribution of P
# and PI is the diagonal matrix of π

# DONE !!!!!!!!!!!!!!!!
hitting_time_digraph <- function(transition)
{
  Ones <- rep(1,dim(transition)[1])
  markovB<-new("markovchain", states=dimnames(transition)[[1]], transitionMatrix=transition)
  pi_<-steadyStates(markovB)[1,]
  PI <- diag(pi_)
  L <- PI - PI%*%transition
  Z <- ginv(L + pi_%*%t(pi_))
  H <- Ones%*%t(diag(Z))-Z
H
}


genetDistStepping <- function(migration,popSize,mutation_rate)
{
  commute_time <- commute_time_undigraph(migration)
  #genetic_dist = commute_time / (8* popSize)
  #genetic_dist = commute_time / (8* (sum(popSize)/(dim(popSize)[1]*dim(popSize)[2])))
  genetDist = (commute_time / 4 + 2*sum(valuesA(popSize)) )* mutation_rate
  genetDist
}


# Calcul of genetic distance from resistance undirected graph

# DONE!!!!!!!!!!!!!!!!
genetDistUndigraph <- function(transition,popSize,mutation_rate)
{
  commute_time <- commute_time_undigraph(transition)
  #genetic_dist = commute_time / (8* popSize)
  #genetic_dist = commute_time / (8* (sum(popSize)/(dim(popSize)[1]*dim(popSize)[2])))
  genetDist =  (commute_time/4+2*sum(valuesA(popSize))) * mutation_rate
  genetDist
}

genetDistUndigraphForNLM <- function(parameters)
{
  migration <- migrationMatrix(rasterStack = environmentalData,shapeDisp = prior$dispersion$model,pDisp = parameters["dispersion"])
  rasK=environmentalData;rasK[cellNumA(rasK)]= as.matrix(ReactNorm(valuesA(environmentalData),p=parameters["K"],shapes=prior$K$model)[,"Y"]);rasK=rasK[[1]]
  transition <- transitionMatrixBackward(r = parameters["R"],K = valuesA(rasK), migration)
  genetDist_Undigraph <- genetDistUndigraph(transition,popSize=rasK,mutation_rate = parameters["mutation_rate"])
  sum((as.dist(genetDist)-as.dist(genetDist_Undigraph))^2+
        (diag(genetDist)-diag(genetist_Undigraph))^2)/
    (nrow(genetDist)*(nrow(genetDist)+1)/2)
}

linearizedFstUndigraph <- function(transition, popSize)
{
  commute_time <- commute_time_undigraph(transition)
  #genetic_dist = commute_time / (8* popSize)
  #genetic_dist = commute_time / (8* (sum(popSize)/(dim(popSize)[1]*dim(popSize)[2])))
  linearizedFst = commute_time / (16*sum(valuesA(popSize))*ncellA(popSize))
  linearizedFst
}

Collisionijk <- function(Hitting_mat)# expected first collision time of i and j in k
{  
  Tijk=array(NA,dim=c(dim(Hitting_mat)[1],dim(Hitting_mat)[2],dim(Hitting_mat)[1]))
  for (k in 1:dim(Hitting_mat)[1]){
    for (i in 1:dim(Hitting_mat)[1]){
      for (j in 1:dim(Hitting_mat)[2]){
        Tijk[i,j,k] <- max(Hitting_mat[i,k],Hitting_mat[j,k]) 
      }
    }
  }
Tijk
}

# Calcul of genetic distance from digraph
linearizedFstDigraph <- function(transition, popSize)#popSize is raster class
{
  H <- hitting_time_digraph(transition)
  Tijk <- Collisionijk(H)
  Collision <- H
  for (i in 1:dim(Tijk)[1]){
    for (j in 1:dim(Tijk)[2]){
      Collision[i,j] <- Tijk[i,j,]%*%t(diag(transition)/t(valuesA(popSize))/sum(diag(transition)/t(valuesA(popSize))))
    }
  }
Collision
}

# Calcul of genetic distance from digraph
linearizedFstDigraph <- function(transition, popSize)#popSize is raster class
{
  H <- hitting_time_digraph(transition)
  dim2 <- dim(H);dim2[[3]]=2
  H2 <- array(c(H,t(H)),dim=dim2)
  MaxH <- apply(H2,c(1,2),max)
  #H2[,,1] <- MinH; H2[,,2] = (H+t(H))/2
  #MinH2 <- apply(H2,c(1,2),min)
  genetic_dist = MaxH / (8*sum(valuesA(popSize))*ncellA(popSize))
  genetic_dist
}

# DONE!!!!!!!!!!!!
genetDistDigraph <- function(transition,popSize,mutation_rate,method="Goldstein95")
{
  H <- hitting_time_digraph(transition)
  dim2 <- dim(H);dim2[[3]]=2
  H2 <- array(c(H,t(H)),dim=dim2)
  MinH <- apply(H2,c(1,2),min)
  #H2[,,1] <- MinH; H2[,,2] = (H+t(H))/2
  #MinH2 <- apply(H2,c(1,2),min)
  genetic_dist = (MinH/2 + 2*sum(valuesA(popSize)))* mutation_rate
  genetic_dist
}

# Calcul of genetic distance from digraph
linearizedFstDigraph1 <- function(transition, popSize)#popSize is raster class
{
  H <- hitting_time_digraph(transition)
  Coalij <- H
  Coal = array (NA,dim=list(dim(H)[1],dim(H)[2],dim(H)[2]))
  for (i in 1:dim(Coal)[1]){
    for (j in 1:dim(Coal)[2]){
      for (k in 1:dim(Coal)[3]){
        Coal[i,j,k] <- max(H[i,k],H[j,k])
      }
      Coalij[i,j] <- mean(Coal[i,j,])
    }
  }
  linearizedFst = Coalij / (2*sum(valuesA(popSize))*ncellA(popSize))
  linearizedFst
}

comuteTimeDigraph <- function(transition, popSize)#popSize is raster class
{
  H <- hitting_time_digraph(transition)
  commuteTime = (H+t(H)) 
}

#
# distribution of coalescence time among pair of genes depending on their deme of sampling at t=0
# arguments : max_time_interval : maximum time interval to jump over generations
#                                 in the probability calculation
#
#
# value: list 
#        -first : array of coalescence times for deme pairs
#        -second : matrix of expected coalescence times (mean of the distribution)
#

coalescence_prob_time_distribution_matrix <- function(transition,max_time_interval=4,rasK,threshold=1E-6)
{
  # calculates the probability distribution of coalescence among pairs of
  # genes in different cells (cell 1 as line, cell 2 as column)
  # knowing "transition" = gene backward transition probability among demes
  #         "time"
  #         "max_time_interval" = maximum number of generations to group 
  #                           the computation of coalescence
  #         "rasK" = population sizes
  # uses function "matrix.power" of package matrixcalc
  # value:
  # a list of matrix of coalescence probability among genes in the deme graph 
  # at different times in the past starting in the following sequence
  # 1,2,4,16,32,...,2^floor(log2(max_time_grouping)),2*2^floor(log2(max_time_grouping))
  # 3*2^floor(log2(max_time_grouping)) ... n*2^floor(log2(max_time_grouping))
  # untill probability of coalescence is lower than a threshold and reducing everywhere in 
  # the matrix
  time_points <- NA # this vector will contain all the time points where probabilities
                    # will be calculated
  occurence_prob <- array(NA, dim=c(dim(transition),1))# occurence probabilities among demes of individuals at the generation
                          # given in the time_points values at the same position
                          # in time_point vector as in the list
  occurence_and_coalescence <- array(NA, dim=c(dim(transition),1))# occurence and coalescence probabilities among demes of 
                                     # second  individual at the time considered in the 
                                     # time_point vector
  popSizes_receiving <- matrix(valuesA(rasK),nrow=dim(transition)[1],ncol=dim(transition)[2],byrow=TRUE)
  coal_column <- 1/(2*popSizes_receiving)
  coal_column[coal_column>1]<-1
  coal_column_tiMax <- 1-(1-coal_column)^max_time_interval
  not_coalesced_prob <- array(NA, dim=c(dim(transition),1))
  coalescence_prob <- array(NA, dim=c(dim(transition),1))#array(NA,dim=)
  time_points <- 1
  occurence_prob[,,1] <-  diag(1,dim(transition)[1])
  occurence_and_coalescence[,,1] <- occurence_prob[,,1] / (2*popSizes_receiving)
  coalescence_prob[,,1] <- occurence_prob[,,1] %*% t(occurence_and_coalescence[,,1])
  coalescence_prob[,,1][coalescence_prob[,,1]>1]=1
  not_coalesced_prob[,,1] <- 1-coalescence_prob[,,1]
  condition=TRUE
  i=1
  maxtransition <- matrix.power(transition,max_time_interval)
  while (condition){
    i=i+1
    time_interval <- time_points[i-1]
    if (time_interval>max_time_interval) {
      time_interval <- max_time_interval;transition_powered = maxtransition
      coal_column_ti <- coal_column_tiMax
    } else {
      transition_powered = matrix.power(transition,time_interval)
      coal_column_ti <- 1-(1-coal_column)^time_interval
    }
    time_points[i] <- time_points[i-1]+time_interval
    occurence_prob <- abind(occurence_prob,occurence_prob[,,i-1] %*% transition_powered,along=3)
    occurence_and_coalescence <- abind(occurence_and_coalescence,occurence_prob[,,i] *coal_column_ti,along=3)
    coalescence_prob <- abind(coalescence_prob,occurence_prob[,,i] %*% t(occurence_and_coalescence[,,i]),along=3)
    coalescence_prob[,,i][coalescence_prob[,,i]>1]=1
    not_coalesced_prob <- abind(not_coalesced_prob,not_coalesced_prob[,,i-1]*(1-coalescence_prob[,,i]),along=3)
    coalescence_prob[,,i] <- coalescence_prob[,,i]*not_coalesced_prob[,,i-1]
    #condition = any(coalescence_prob[[i]]>=coalescence_prob[[i-1]])|any(coalescence_prob[[i]]>threshold)
    condition = !all(not_coalesced_prob[,,i]<threshold)
  }
  dimnames(coalescence_prob) <- list(1:dim(transition)[1],1:dim(transition)[2],time_points)
  expected_coalescence_times <- matrix(NA,nrow=dim(transition)[1],ncol=dim(transition)[2])
  for (i in 1:dim(transition)[1]){
    for (j in 1:dim(transition)[2]){
      expected_coalescence_times[i,j] <- coalescence_prob[i,j,]%*%time_points
    }
  }
  #expected_coalescence_times<-array(unlist(coalescence_prob),dim=c(dim(transition),length(coalescence_prob)))
  #expected_coalescence_times<-unlist(coalescence_prob,dim=c(dim(matrix),length(coalescence_prob))
  list(coalescent_prob=coalescence_prob,exp_times=expected_coalescence_times)
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
#

CreateGenetArray <- function(rasK, nb_locus, initial_locus_value,Option="sample_1_col_diploid",nind="Ne")
{
  if (nind%in%c("Ne","Ne2","Ne10")) nind <- switch(nind,
                                                   Ne = floor(sum(valuesA(rasK))),
                                                   Ne2= 2*floor(sum(valuesA(rasK))),
                                                   Ne10= 10*floor(sum(valuesA(rasK)))
  )
  #Get coords for each cell
  coords = xyFromCellA(rasK)
  repet = switch(Option,
                 sample_1col_diploid= sort(rep(sample(rep(1:ncellA(rasK),round(valuesA(rasK))),nind,replace=TRUE),2)), 
                 sample_2col_diploid= sample(rep(1:length(rasK),round(valuesA(rasK))),nind,replace=TRUE), 
                 sample_haploid= sample(rep(1:ncellA(rasK),round(valuesA(rasK))),nind,replace=TRUE),
                 full_1col_diploid= rep(1:length(rasK),round(valuesA(rasK))*2), 
                 full_2col_diploid= rep(1:length(rasK),round(valuesA(rasK))), 
                 full_haploid= rep(1:length(rasK),round(valuesA(rasK))),
                 homogenous_1col= rep(1:length(rasK),each=nind/length(rasK)*2),
                 homoNonEmtpyCells1col= rep((1:length(rasK))[which(round(valuesA(rasK))>1)],each=nind/length(rasK)*2)
  )
  geneticData <- as.data.frame(coords[repet,]) ; geneticData[,"Cell_numbers"] <- repet
  genes = switch(Option,
                 sample_1col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                 sample_2col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                 sample_haploid =      as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                 full_1col_diploid =   as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                 full_2col_diploid =   as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                 full_haploid =        as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                 homogenous_1col =     as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                 homogenous_2col =     as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus)*2)),
                 homoNonEmtpyCells1col=as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus)))
  )
  colnames(genes) = switch(Option,
                           sample_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                           sample_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                           sample_haploid = paste("Locus",1:nb_locus, sep=""),
                           full_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                           full_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                           full_haploid = paste("Locus",1:nb_locus, sep=""),                 
                           homogenous_1col = paste("Locus",1:nb_locus, sep=""),                 
                           homoNonEmtpyCells1col= paste("Locus",1:nb_locus, sep=""),                 
                           homogenous_2col =  paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep="")
  )
  geneticData <- cbind(geneticData,genes) # add locus to geneticData 
  geneticData
}



#Plot genetic data in environmental data observed 
#Genetic data is turn into a Spatial Pixel Data Frame 
#Mettre des couleurs en fonction du nombre d'individu
plotGeneticData = function(geneticData, EnvironmentalDataObserved)
{
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
repnDispMutFunction <- function(geneticData, dimGeneticData, mutationRate, transitionmatrice)
{
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
repnDispMutFunction <- function(geneticData, dimGeneticData, mutationRate, transitionmatrice)
{
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
Qwithin_pair <- function(gDat)
{
  gDat <- gDat[,grep("ocus",colnames(gDat))]
  if ((dim(gDat)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat)))) { #gDat is 1_col_diploid
    gDat <- OneCol2TwoCols(gDat)}
  matrix_pair = gDat[,grep(".2", colnames(gDat), fixed = T)]
  matrix_impair = gDat[,grep(".1", colnames(gDat), fixed = T)]
  Qw <- (matrix_pair == matrix_impair)
  Qw = rowMeans(Qw) # vector of probability of Qw for each individual
  Qw = matrix(Qw, ncol = dim(gDat)[1], nrow = dim(gDat)[1])
  Qw = (Qw+t(Qw))/2
Qw
}

Qwithin_pair_2 <- function(gDat)
{
  gDatG <- gDat[,grep("ocus",colnames(gDat))]
  gDat3D <- OneCol23Dims(gDatG)
  (t(array(diag(id_mat_prob(gDat3D[,,1],gDat3D[,,2])),dim=c(dim(gDat3D)[1],dim(gDat3D)[1])))+array(diag(id_mat_prob(gDat3D[,,1],gDat3D[,,2])),dim=c(dim(gDat3D)[1],dim(gDat3D)[1])))/2
}

#Function that calculate probability of identity of genes intra individual (at population level)
Qwithin_pop <- function(gDat)
{
  gDat <- gDat[,grep("ocus",colnames(gDat))]
  if ((dim(gDat)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat)))) { #gDat is 1_col_diploid
  gDat <- OneCol2TwoCols(gDat)}
  matrix_pair = gDat[,grep(".2", colnames(gDat), fixed = T)]
  matrix_impair = gDat[,grep(".1", colnames(gDat), fixed = T)]
  Qw <- (matrix_pair == matrix_impair)
  Qw = mean(Qw)
  Qw
}

OneCol2TwoCols <- function(gDat)
{
  NonGeneticCols <- gDat[,-grep("ocus",colnames(gDat))]
  gDat <- gDat[,grep("ocus",colnames(gDat))]
  if (dim(gDat)[1] %% 2 != 0) stop("odd number of rows in suposedly one col genetic table")
  pair <- (1:(dim(gDat)[1]/2))*2;impair=pair-1
  first <- gDat[impair,]
  colnames(first) <- paste(colnames(first),"1",sep=".")
  second <- gDat[pair,]
  colnames(second) <- paste(colnames(second),"2",sep=".")
  Twocols<-cbind(first,second)
  neworder <- rep(c(0,dim(gDat)[2]),dim(gDat)[2])+rep(1:dim(gDat)[2],each=2)
  Twocols <- Twocols[,neworder]
  cbind(NonGeneticCols[impair,],Twocols)
}

TwoCols2OneCol_ <- function(gDat)
{
  notlocusonly <- (length((grep("ocus",colnames(gDat))))<dim(gDat)[2]) 
  if (notlocusonly) {
    coords <- gDat[rep(1:dim(gDat)[2],each=2),-grep("ocus",colnames(gDat))]
    coords <- data.frame(x=rep(gDat[,"x"],each=2),y=rep(gDat[,"y"],each=2),Cell_numbers=rep(gDat[,"Cell_numbers"],each=2))
    rownames(coords)<-paste(paste(rep(rownames(gDat),each=2),each=1:2,sep="."))
  } 
  gDat <- gDat[,grep("ocus",colnames(gDat))]
  pair <- (1:(dim(gDat)[2]/2))*2;impair=pair-1
  first <- gDat[,impair];colnames(first) <- unlist(strsplit(colnames(first),"\\."))[impair]
  second <- gDat[,pair];colnames(second)<-colnames(first)
  Onecol<-rbind(first,second)
  neworder <- rep(c(0,dim(gDat)[1]),dim(gDat)[1])+rep(1:dim(gDat)[1],each=2)
  Onecol <- Onecol[neworder,]
  row.names(Onecol) <- paste(paste(rep(rownames(gDat),each=2),each=1:2,sep="."))
  if (notlocusonly) (cbind(coords,Onecol)) else Onecol
}

TwoCols2OneCol <- function(tip_genotype)
{
  nonGenotypeData <- tip_genotype[,!grepl("ocus",colnames(tip_genotype))]
  
  matrix_pair = tip_genotype[,grep("\\.2", colnames(tip_genotype))]
  colnames(matrix_pair) <- grep("ocus",unlist(strsplit(colnames(matrix_pair),"\\.")),value=TRUE)
  matrix_pair <- cbind(nonGenotypeData,matrix_pair)
  rownames(matrix_pair)<-paste(rownames(matrix_pair),".1",sep="")
  matrix_impair = tip_genotype[,grep("\\.1", colnames(tip_genotype))]
  colnames(matrix_impair) <- grep("ocus",unlist(strsplit(colnames(matrix_impair),"\\.")),value=TRUE)  
  matrix_impair <- cbind(nonGenotypeData,matrix_impair)
  rownames(matrix_impair)<-paste(rownames(matrix_impair),".2",sep="")  
  tip_genotype <- rbind(matrix_pair,matrix_impair)[rep(1:nrow(matrix_pair),each=2)+rep(c(nrow(matrix_pair),0),nrow(matrix_pair)),]
  tip_genotype
}

OneCol23Dims <- function(gDat)
{
  first <- gDat[2*(1:(dim(gDat)[1]/2))-1,grep("ocus",colnames(gDat))]  
  second <- gDat[2*(1:(dim(gDat)[1]/2)),grep("ocus",colnames(gDat))]  
  dimN <- dimnames(first);dimN[[3]]<-c("al.1","al.2")
  array(c(unlist(first),unlist(second)),dim=c(dim(first),2),dimnames=dimN)
}

id_mat_prob <- function(A,B)
{
  id_mat_prob <- array(NA,dim=c(dim(A)[1],dim(A)[1]),dimnames=list(unlist(dimnames(A)[1]),unlist(dimnames(A)[1])))
  for (i in 1:dim(A)[1]){
    for (j in 1:dim(A)[1]){
      id_mat_prob[i,j] <- mean(A[i,]==B[j,])
    }
  }
  id_mat_prob
}

Qbetween_2 <- function(gDat)
{
  gDatG <- gDat[,grep("ocus",colnames(gDat))]
  gDat3D <- OneCol23Dims(gDatG)
  (id_mat_prob(gDat3D[,,1],gDat3D[,,1])
  +id_mat_prob(gDat3D[,,2],gDat3D[,,2])
  +id_mat_prob(gDat3D[,,1],gDat3D[,,2])
  +id_mat_prob(gDat3D[,,2],gDat3D[,,1]))/4
}

# Fonction that calculates probability of identity of genes inter individual (between two individuals)
Qbetween <- function(gDat)
{
  colNames  <- colnames(gDat)
  gDat <- gDat[,grep("ocus",colnames(gDat))]
  if ((dim(gDat)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat)))) { #gDat is 1_col_diploid
    gDat <- OneCol2TwoCols(gDat)}
  Qb = matrix(ncol = dim(gDat)[1], nrow = dim(gDat)[1]) #initialization of Qb as a matrix
  # A = genetic data with loci only
  A=as.matrix(gDat[,grep("ocus",colnames(gDat),fixed=T)])
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

SSb <- function(gDat)
{
  1+Qwithin_pair(gDat)-2*Qbetween(gDat)
}

SSw <- function(gDat)
{
  2*(1-Qwithin_pop(gDat)) 
}

a_matrix <- function(gDat)
{
  SSb(gDat)/SSw(gDat)-1/2
}

# TEST: pour un nombre de generation donné, on teste la stabilité du a value
#Function test of stabilisation for a value
test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice)
{
## ref: Rousset et al. J Evol Biol,13  (2000) 58-62.
  vecteur_a_value <-c(0)
  for(i in 1: nb_generations){
    print(i)
    geneticData <- repnDispMutFunction(geneticData, dimGeneticData, mutationRate, transitionmatrice) 
    matrixQb = (1-Qwithin_pair(geneticData)+2*(Qwithin_pair(geneticData)-Qbetween(geneticData)))
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

test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice)
{
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


# this fstat function is not finished and not used
fstat = function(geneticData)
{
  Genotypes = geneticData[,grep("Locus", colnames(geneticData), fixed = T)]
  form <- as.formulae
  Pops = geneticData[,"Cell_numbers"]
  MeanPop = t(matrix((colSums(Genotypes)/dimGeneticData[1]),ncol=dimGeneticData[1],nrow=dimGeneticData[2]))
  VarInd = matrix(Genotypes^2 - MeanTot^2,ncol=dimGeneticData[2])
  VarInterPop = var(MeanPop)
  VarIntraPop = colSums(VarInd)/dimGeneticData[1]
  VarTot = VarInd
}

# this mmute function is not used
mmute <- function(cells=c(1,2),transitionmatrice)
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
# this ocur function is not used
ocur <- function(cells=c(1,2),transitionmatrice)
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

#
# transitionMatrice 
# function that calculates matrixes of transition, forward and backwark, as well as carrying capacity
# depending on environmentall data and r, K niche model parameters and shapes
#
transitionMatrice <- function(rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp)
{
  K = ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"]
  r = ReactNorm(valuesA(rasterStack),pR,shapesR)[,"Y"]
  migrationM <- migrationMatrix(rasterStack,shapeDisp, pDisp)
  transitionmatrice = transitionMatrixBackward(r, K, migration= migrationM);
  transition_forward = transitionMatrixForward(r, K, migration= migrationM)
  list(backw=transitionmatrice,forw=transition_forward,K=K)  
}

# coalescent
# function that simulates a coalescent in a lansdcape characterized by an environmental variable raster
# stack in aspecies with given niche function relating the carrying capacity and growth rate 
# with the environemental variable 
# pK : parameters of K / environmental variables
# pR : parameters of r / environmental variables
# shapesK : shapes of niche funciotn (reaction norm)
# pDisp : parameters of dispersion 
# rasterStack : environmental variables raster stack
# genetic data table : with coordinates
# aggre_gener : number of generations to aggregate in the simulation steps

simul_coalescent <- function(transitionList,geneticData)
{
  prob_forward=NA
  N <- round(transitionList$K);N[N==0]<-1
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
  for (cell in 1:ncellA(rasterStack))#cell=1)
  {
    nodes_remaining_by_cell[[cell]] <- which(cell_number_of_nodes==cell)
  }
  while (length(unlist(nodes_remaining_by_cell))>1) 
  {
    # migration
    # we localize the parents in the landscape by sampling in the backward transition matrix
    for (node in 1:length(parent_cell_number_of_nodes))#gene=1;node=1# parent_cell_number_of_nodes
    {
      parent_cell_number_of_nodes[node] = sample(ncellA(rasterStack),size=1,prob=c(transitionList$backw[cell_number_of_nodes[node],]))
    }
    # once we know the parent cell numbers, we calculate the forward dispersion probability of the event
    prob_forward[time] = sum(log(transitionList$forw[parent_cell_number_of_nodes,cell_number_of_nodes]))
    # coalescence
    time=time+1; if (round(time/10)*10==time) {print(time)}
    # we now perform coalescence within each cell of the landscape for the parents
    for (cell in 1:ncellA(rasterStack))#cell=1;cell=2;cell=3;cell=4;cell=5;cell=26;cell=10
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
  list(coalescent=coalescent,prob_forward=sum(prob_forward))
}

#
# simul_coal_genet
# simulates coalescent depending on environmental layers and r/K niche funciton parameters 
# and simulates genetic evolution along the coalescent
# 

simul_coal_genet <- function(geneticData,rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp,mutation_rate=1E-2,initial_locus_value=200,mutation_model="stepwise",stepvalue=2,locusnames=NA)
{
  if (is.na(locusnames)) {
    locusnames <- grep("ocus",colnames(geneticData),value=TRUE)
    # if locusnames is of the type loc1.1 loc1.2 loc2.1 loc.2.2 ...: already 2col data
    if (sum(grep("\\.1",locusnames,value=FALSE)|
             grep("\\.2",locusnames,value=FALSE))==length(locusnames)) {
      allelenames<-locusnames} else {
        allelenames <- paste(rep(locusnames,each=2),".",1:2,sep="")
      }                           
  }
  tip_genotype <- geneticData[,1:3]
  coalescent <- list()
  transitionMatrixList <- transitionMatrice(rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp)
  for (allele in allelenames){
    coalescent[[allele]] <- simul_coalescent(transitionMatrixList,geneticData)
    coalescent_mut=add_br_length_and_mutation(coalescent[[allele]]$coalescent,mutation_rate,initial_locus_value,allele)
    coaltable <- coalist_2_coaltable(coalescent_mut,allele,initial_locus_value)
    genetic_values=genetics_of_coaltable(coaltable,initial_locus_value,mutation_model,stepvalue)
    rownames(genetic_values) <- genetic_values$coalescing
    tip_genotype[,allele]<-genetic_values[as.character(1:dim(geneticData)[1]),"genotype"]
  }
  list(coalescent=coalescent,mutation_rate=mutation_rate,tip_genotype=tip_genotype,
       transitionAndK=transitionMatrice(rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp))
  # forward_log_prob is the average per generation of the log probability of the forward movements of the genes in the landscape
}

#
# coalist_2_coaltable function converts a coalescent list to a coalecent table format
#

coalist_2_coaltable <- function(coalist,locusnames,initial_locus_value)
{
  coaldf <- data.frame(Reduce(rbind,coalist))
  coaltable <- coaldf[rep(1:dim(coaldf)[1],unlist(lapply(coaldf$coalescing, length))),]
#  coaltable[,c("coalescing","br_length","mutations")] <- unlist(coaldf[,c("coalescing","br_length","mutations")])
  coaltable[,c("coalescing","br_length",locusnames)] <- unlist(coaldf[,c("coalescing","br_length",locusnames)])
  coaltable$genotype <- NA
  coaltable[dim(coaltable)[1]+1,"genotype"] <- initial_locus_value
  coaltable[dim(coaltable)[1],"coalescing"] <- dim(coaltable)[1]
  coaltable$mut <- coaltable[,grep("Locus",colnames(coaltable),value=TRUE)]
  coaltable <- coaltable[,-which(colnames(coaltable)%in%grep("Locus",colnames(coaltable),value=TRUE))]
  #rownames(coaltable) <- coaltable$coalescing
}

#
# genetics_of_coaltable function
# adds genetic values to a coalescent table containing mutation number per branch
# knowing initial genetic value of the ancastor and mutation model

genetics_of_coaltable <- function(coaltable,initial_locus_value,mutation_model="stepwise",stepvalue=2)
{
 switch(mutation_model,
        stepwise = stepwise(coaltable,initial_locus_value,stepvalue,locusnames),
        my_ass = "rien"
        ) 
 #stepwise(coaltable,initial_locus_value,stepvalue, locusnames)
}

#
# stepwise : simulates stepwise mutation process along a coalescent table  
#            using the number of mutation per branch information ("mut" column of the 
#            coaltable)
#
stepwise <- function(coaltable,initial_locus_value,stepvalue,locusnames)
{
  #  coaltable$genetic_value=NA
  # we calculate the oritattion of the mutations in the different branches using binomial rules
  #coaltable[dim(coaltable)[1]+1,] <- c(NA,max(unlist(coaltable$new_node)),NA,NA,NA,initial_locus_value,NA)
  #coaltable$directional = 2*rbinom(prod(dim(coaltable[,"mut"])),unlist(coaltable[,"br_length"]*coaltable[,"mut"]),.5)-coaltable[,"br_length"]*coaltable[,"mut"]
  coaltable$dir = 2*rbinom(n=prod(dim(as.table(coaltable[,"mut"]))),unlist(coaltable[,"mut"]),.5)-coaltable[,"mut"]
  coaltable$genotype <- is.na(coaltable$mut)
  coaltable[dim(coaltable)[1],"genotype"] <- initial_locus_value
  #  for (locus in locusnames)
  #  {  
  for(branch in rev(rownames(coaltable)[-dim(coaltable)[1]]))
  {
    coaltable[branch,"genotype"] <- as.matrix(coaltable[branch,"dir"]*stepvalue + coaltable[which(coaltable$coalescing==coaltable[branch,"new_node"]),"genotype"])
  }
  #  }
  coaltable
}

#
#
#
stepwise <- function(coaltable,initial_locus_value,stepvalue,locusnames)
{
  #  coaltable$genetic_value=NA
  # we calculate the oritattion of the mutations in the different branches using binomial rules
  #coaltable[dim(coaltable)[1]+1,] <- c(NA,max(unlist(coaltable$new_node)),NA,NA,NA,initial_locus_value,NA)
  #coaltable$directional = 2*rbinom(prod(dim(coaltable[,"mut"])),unlist(coaltable[,"br_length"]*coaltable[,"mut"]),.5)-coaltable[,"br_length"]*coaltable[,"mut"]
  coaltable$dir = 2*rbinom(n=prod(dim(as.table(coaltable[,"mut"]))),unlist(coaltable[,"mut"]),.5)-coaltable[,"mut"]
  coaltable$genotype <- is.na(coaltable$mut)
  coaltable[dim(coaltable)[1],"genotype"] <- initial_locus_value
  #  for (locus in locusnames)
  #  {  
  for(branch in rev(rownames(coaltable)[-dim(coaltable)[1]]))
  {
    coaltable[branch,"genotype"] <- as.matrix(coaltable[branch,"dir"]*stepvalue + coaltable[which(coaltable$coalescing==coaltable[branch,"new_node"]),"genotype"])
  }
  #  }
  coaltable
}
#
# add br_length and mutation to coalescent list
#
#

add_br_length_and_mutation <- function(coalescent,mutation_rate,initial_locus_value,allelenames)
{
  tips = NULL
  internals = NULL
  nodes = NULL
  times = NULL
  for (i in 1:length(coalescent))#i=1;i=2
  {
    nodes = append(nodes,c(coalescent[[i]]$coalescing,coalescent[[i]]$new_node))
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
      for (allele in allelenames)
      {
        coalescent[[i]][[allele]] <- rpois(rep(1,length(coalescent[[i]]$br_length)),coalescent[[i]]$br_length*mutation_rate)
      }
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
# goldstein deltamu^2 genetic distance
#

deltaMu <- function(tip_genotype,stepvalue=2)
{
  g <- aggregate(tip_genotype,by=list(tip_genotype$Cell_numbers),FUN="mean",na.action="na.omit")
  Cells <-g$Cell_numbers  
  g <- g[,grep("ocus",colnames(tip_genotype))]
  gd <- as.matrix(dist(g/stepvalue)^2/dim(g)[2])
  dimnames(gd) <- list(Cells,Cells)
  gd
}

#
# Shared allele distance
#

SharedAlleleDistance <- function(tip_genotype)
{
  if (all(grep("\\.2",colnames(tip_genotype))==(grep("\\.1",colnames(tip_genotype))+1))) {
    tip_genotype <- TwoCols2OneCol(tip_genotype)
    } 
  Cells <- levels(as.factor(tip_genotype$Cell_numbers))
  shAll_mat <- matrix(NA, nrow=length(Cells),ncol=length(Cells))
  for (cell1 in Cells)
    for(cell2 in Cells)
    {
      alleles1 <- tip_genotype[tip_genotype$Cell_numbers==cell1,]
      alleles2 <- tip_genotype[tip_genotype$Cell_numbers==cell2,]
    }
  g <- tip_genotype[,grep("ocus",colnames(tip_genotype))]
  g1 <- g[,1:ncol(g)/2]
  g2 <- g[,1+1:ncol(g)/2]
  #g[] <- 
  #
  # ! not finished
}

#
# genetic distance among nodes in a genotype table
#

genetDist <- function(tip_genotype,method="deltaMu",stepvalue=2,byCell=TRUE)
{
  genetDist =switch(method,
                    deltaMu=deltaMu(tip_genotype,stepvalue=stepvalue))
genetDist
}

get_nj_tree <- function(tip_genotype,mutation_model,step_value)
{
  if (mutation_model=="stepwise") distmat <- deltaMu(tip_genotype,stepvalue=step_value)
  nj(as.matrix(distmat))
}

get_ultrametric_nj_tree <- function(tip_genotype,mutation_model,step_value)
{
  #
  nj_tree <- get_nj_tree(tip_genotype,mutation_model,step_value)
  nj_tree$edge.length[nj_tree$edge.length<0]<-0 # negative branch length set to zero
  chronos(nj_tree)
}

#
# Coalescent probabilities
#



#
# tmra of coalescent : tiem to most recent comon ancestor
#

tmra <- function(coalescent_simulated)
{
  coalescent_simulated$coalescent[[length(coalescent_simulated$coalescent)]]$time
}
#
# plot_coalescent plots a coalescent simulation
# argument: output of alescent()

plot_genealogy <- function(genealogy,file=NA)
{
  if (!is.na(file)) pdf(file)
  par(mfrow=c(2,1),oma=c(0,0,0,4),xpd=TRUE)
#  tipcells <- geneticData$Cell_numbers[as.numeric(coalescent_2_phylog(genealogy)$tip.label)]
#  tipcols = rainbow(ncellA(rasK))[tipcells]
  plot(coalescent_2_phylog(genealogy),direction="downward")
#  legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(legend_right_move,0))
#  text(x=7,y=0,paste("tmra = ",tmra(coalescent_simulated)))
#  njtree<-get_nj_tree(coalescent_simulated$tip_genotype,mutation_model="stepwise",step_value=2)
#  njtree$edge.length[njtree$edge.length<0]<-0 # negative branch length set to zero
#  njtree$edge.length<-njtree$edge.length/10^floor(log10(max(njtree$edge.length))) # scale tree to plot
#  plot(njtree,direction="downward")
#  if (with_landscape) {plot(rasK)}
  if (!is.na(file)) dev.off()
  
}

#plot_coalescent old one
plot_coalescent <- function(coalescent_simulated,with_landscape=FALSE,rasK=NULL,legend_right_move=-.2,file=NA)
{
  coalescent <- coalescent_simulated$coalescent
  if (!is.na(file)) pdf(file)
  if (with_landscape) {par(mfrow=c(3,1),oma=c(0,0,0,4),xpd=TRUE)}else{par(mfrow=c(2,1),oma=c(0,0,0,4),xpd=TRUE)}
  tipcells <- geneticData$Cell_numbers[as.numeric(coalescent_2_phylog(coalescent)$tip.label)]
  tipcols = rainbow(ncellA(rasK))[tipcells]
  plot(coalescent_2_phylog(coalescent),direction="downward",tip.color=tipcols)
  legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(legend_right_move,0))
  text(x=7,y=0,paste("tmra = ",tmra(coalescent_simulated)))
  njtree<-get_nj_tree(coalescent_simulated$tip_genotype,mutation_model="stepwise",step_value=2)
  njtree$edge.length[njtree$edge.length<0]<-0 # negative branch length set to zero
  njtree$edge.length<-njtree$edge.length/10^floor(log10(max(njtree$edge.length))) # scale tree to plot
  plot(njtree,direction="downward")
  if (with_landscape) {plot(rasK)}
  if (!is.na(file)) dev.off()
}

#
#plot_coalescent new one
#

plot_coalescent <- function(coalescent_simulated,with_landscape=FALSE,rasK=NULL,legend_right_move=-.2,file=NA)
{
  coalescent <- coalescent_simulated$coalescent
  gDat <- coalescent_simulated$tip_genotypes
  if (!is.na(file)) pdf(file)
  if (with_landscape) {par(mfrow=c(3,1),oma=c(0,0,0,4),xpd=TRUE)}else{par(mfrow=c(2,1),oma=c(0,0,0,4),xpd=TRUE)}
  tipcells <- gDat$Cell_numbers[as.numeric(coalescent_2_phylog(coalescent)$tip.label)]
  tipcols = rainbow(ncellA(rasK))[tipcells]
  plot(coalescent_2_phylog(coalescent),direction="downward",tip.color=tipcols)
  legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(legend_right_move,0))
  text(x=7,y=0,paste("tmra = ",tmra(coalescent)))
  njtree<-get_nj_tree(coalescent,mutation_model)
  njtree$edge.length[njtree$edge.length<0]<-0 # negative branch length set to zero
  njtree$edge.length<-njtree$edge.length/10^floor(log10(max(njtree$edge.length))) # scale tree to plot
  plot(njtree,direction="downward")
  if (with_landscape) {plot(rasK)}
  if (!is.na(file)) dev.off()
}
#
# geneticDataSimulList : a list of simulations, with sublist geneticData and sublist log_lik_forward 
#

summary_stat <- function(geneticDataObs,geneticDataSimulList,log_lik_simul_list)
{
  #1) We calculate a matrix of observed genetic distance between individuals (DSAI = sum_j (Shared alleles) / (number of loci))
  #   or reapeat number distance (Goldstein 1995)
  #2) We select PCi representing 99% cumulative variance
  #3) We express simulated data in these axis. First summary stats are the values of each invidiual on these major axes
  #4) Second summary stats are mean number of alleles per individual for observed and simulated data
  #5) Third summary stats are mean number of alleles per population for observed and simulated data
  #
  
}




#
# genetic_simulation
# this function simulates genetic data from a coalescent
#
#


test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice)
{
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

#
# Fst and Fis among cells of a raster given in the gDat genetic table 
#
FstatsRas <- function(gDat,by="cell",all=TRUE,cells=NULL)
{
  Absent <- !cells%in%gDat$Cell_numbers
  if ((all==TRUE)&(any(Absent))){
    AbsentCells <- (cells)[Absent]
    gDat[dim(gDat)[1]+1:length(AbsentCells),"Cell_numbers"] <-AbsentCells
  }
  for (i in levels(as.factor(gDat$Cell_numbers))){
    gDat[gDat$Cell_numbers==as.numeric(i),"newCelN"] <- which(levels (as.factor(gDat$Cell_numbers))==i)
  }
  gDat2 <- gDat[,grep("ocus",colnames(gDat))]
  if ((dim(gDat2)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat2)))) { #gDat is 1_col_diploid
    gDat <- OneCol2TwoCols(gDat)}
  genotypes <- gDat[,grep("ocus",colnames(gDat))]
  cells <- as.integer(levels(as.factor(gDat$Cell_numbers)))
  Fstat(genotypes,npop=length(cells),pop.mbrship=gDat$newCelN,ploidy=2)
}

#
# Genetic_Dist
#

Genetic_Dist <- function(gDat,method="Goldstein")
{
switch(method,
       Goldstein=deltaMu(tip_genotype = gDat,stepvalue = 1)  
)
}

#
# Individual a value (Fst/(1-Fst), Rousset 2000)
#

a_value_ind <- function(gDat)
{
  # gDat is a diploid genotype data frame with Cell_numbers column and 1 or 2 cols per individual
  #
  matrixQb = (1-Qwithin_pair(gDat)+2*(Qwithin_pair(gDat)-Qbetween(gDat)))
  matrixQw = 2*(1-Qwithin_pop(gDat))
  vecteur_a_value = matrixQb/matrixQw-1/2
  vecteur_a_value[is.na(vecteur_a_value)] <-0
  vecteur_a_value
}


# to get a matrix of a-values between all cells of a landscape from 
# simulations of evolution where not all cells necesarily contain individuals
#
#

niche_model <- function(rasterStack,
                        abundanceData,
                        priorsK = list(min=matrix(c(10,200,11,2,0,0,250,1500,1500,2,0,0),
                                    nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),
                                                                c("BIO1","BIO12")))
                                    , max=matrix(c(200,400,300,30,0,0,300,2000,1500,30,0,0),
                                                 nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),
                                                                             c("BIO1","BIO12")))),
                        shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"))
{
  # abundanceData
  # abundanceData$N : number of units sampled
  # abundanceData$k : number of entities observed
  # abundanceData$x : longitude
  # abundanceData$y : latitude
  rasK=rasterStack;valuesA(rasK)= as.matrix(ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"]);rasK=rasK[[1]]
  abundanceData <- cbind(abundancedata,extract(rasK,spatialPoints(abundanceData[,c("x","y")])))
}



a_value_matrix <- function(gDat,rasterStack)
{
  Cell_numbers = gDat$Cell_numbers
  cell_numbers_levels <- levels(as.factor(gDat$Cell_numbers))
  vecteur_a_value <- a_value_ind(gDat)
  gDat <- gDat[,grepl("\\.",colnames(gDat))]
  if ((dim(gDat)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat)))) { #gDat is 1_col_diploid
    cell_numbers_of_ind <- Cell_numbers[(1:(dim(gDat)[1]/2))*2-1]} else {
      cell_numbers_of_ind <- Cell_numbers
    }
  dimnames(vecteur_a_value)<-list(cell_numbers_of_ind,cell_numbers_of_ind)
  agg <- aggregate(vecteur_a_value,by=list(cell_numbers_of_ind),FUN="mean")[,-1]
  agg2 <- t(aggregate(t(agg),by=list(cell_numbers_of_ind),FUN="mean"))[-1,]
  row.names(agg2) <- as.numeric(cell_numbers_levels)
  colnames(agg2) <- as.numeric(cell_numbers_levels)
  agg3 <- matrix(NA, nrow=ncellA(rasterStack),ncol=ncellA(rasterStack))
  agg3[as.numeric(row.names(agg2)),as.numeric(colnames(agg2))] <- agg2
agg3
}

# Function that determine number of dispersal parameters from dispersal shapeDisp used
# Useful for nlm estimation
nbpDisp <- function(shapeDisp)
{
  (switch(shapeDisp,
          fat_tail1 = 2,
          gaussian = 1,
          exponential = 1,
          contiguous = 1,
          island = 1,
          fat_tail2 = 2))
}

# Function that simulate a genetic data with forward simulation, parameters given
#It returns a list of 2 variables:  final genetic data observed and matrix of a-value observed
simulationGenet <- function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, mutationRate, nbLocus, initial_locus_value, shapeDisp, pDisp, nb_generations,ind_per_cell=30)
{
  K <- subset(ReactNorm(valuesA(donneesEnvironmentObs), pK , shapesK),select="Y") # carrying capacity
  r <- subset(ReactNorm(valuesA(donneesEnvironmentObs), pR , shapesR),select="Y") # growth rate
  Rast_K <- donneesEnvironmentObs ; valuesA(Rast_K) <- K
  Rast_r <- donneesEnvironmentObs ; valuesA(Rast_r) <- r
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
  nCell = ncellA(donneesEnvironmentObs)
  Cell_numbers <- cellNumA(donneesEnvironmentObs)
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
expect_linearizedFst_undigraph <- function(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA)
{
  if (is.na(Cell_numbers[1])) Cell_numbers = cellNumA(rasterStack)
  K = ReactNorm(valuesA(rasterStack), pK, shapesK)[,"Y"]
  r = ReactNorm(valuesA(rasterStack), pR, shapesR)[,"Y"]
  rasK=rasterStack;rasK[cellNumA(rasK)]= as.matrix(ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"])
  matrice_migration = migrationMatrix(rasK,shapeDisp, pDisp)
  matrice_transition = transitionMatrixBackward(r,K, matrice_migration)
  linearizedFst = linearizedFstUndigraph(matrice_transition, popSize=rasK)[Cell_numbers,Cell_numbers]
  linearizedFst
}

expect_linearizedFst_digraph <- function(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA)
{
  if (is.na(Cell_numbers[1])) Cell_numbers = cellNumA(rasterStack)
  K = ReactNorm(valuesA(rasterStack), pK, shapesK)[,"Y"]
  r = ReactNorm(valuesA(rasterStack), pR, shapesR)[,"Y"]
  rasK=rasterStack;rasK[cellNumA]= as.matrix(ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"])
  matrice_migration = migrationMatrix(rasK,shapeDisp, pDisp)
  matrice_transition = transitionMatrixBackward(r,K, matrice_migration)
  linearizedFst_att = linearizedFstDigraph(matrice_transition,rasK)[Cell_numbers,Cell_numbers]
  linearizedFst_att
}

ssr <- function(p)
{
  linearizedFst_att = expect_linearizedFst_undigraph(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=finalGenetData$Cell_Number)
  return(mean((linearizedFst_obs - linearizedFst_att)^2))
}

#initial = c(0,1,20,1.5)
#fct_erreur_calc = nlm(f = ssr , p = initial, hessian = FALSE)
#p=initial
#linearizedFst_graph_model = expect_linearizedFst(donneesEnvironmentObs,initial,finalGenetData$Cell_Number,nbpDisp=2,nblayers=1,shapeDisp="fat_tail1")
#ssr(initial)
#parametres_infere = fct_erreur_calc$estimate
#parametres_reels = c(alpha = 0, beta = 2,beta = 1, alphaDisp = 20, betaDisp = 1.5)

#
# validation using forward simulation
#

validation <- function(donneesEnvironmentObs,pK = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))), pR = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))),shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"),shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"), shapeDisp="fat_tail1", pDisp = c(mean=0.32,shape=1.6), file=NULL,mutationRate,nbLocus, initial_locus_value,nb_generations=5000,indpercell=30)
{
  nblayers =dim(donneesEnvironmentObs)[3]
  nCell = ncellA(donneesEnvironmentObs)
  Cell_numbers <- cellNumA(donneesEnvironmentObs)
  K <- ReactNorm(valuesA(donneesEnvironmentObs), pK , shapesK)$Y # recuperation de la carrying capacity
  r <- ReactNorm(valuesA(donneesEnvironmentObs), pR , shapesR)$Y # recuperation du growth rate
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
  transitionmatrice = transitionMatrixBackward(Npop = K, migration= migrationM)
  geneticObs = simulationGenet(donneesEnvironmentObs,alpha, beta, mutationRate,nbLocus, initial_locus_value,shapeDisp,pDisp,nb_generations=5000,indpercell)
  finalGenetData = geneticObs[[1]]
  land_size <- raster(matrix(K[1,],nrow=dim(donneesEnvironmentObs)[1],ncol=dim(donneesEnvironmentObs)[2]))
  a_value_simul = a_value_matrix(geneticObs[[1]],land_size)
  
  a_value_theory_stepping_stone_model = S/(32*dnorm(res(rasK)[1],0,pDisp[1])*sum(valuesA(rasK))/ncellA(rasK))
  distanceMatrix(donneesEnvironmentObs)/(4*0.05)
  a_value_theory_island_model <- matrix(pDisp[1],nrow=nCell,ncol=nCell)-pDisp[1]*diag(nCell)
  a_value_theory_graph_model <- expect_linearizedFst_undigraph(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers)
  if (!is.null(file)){jpeg(file)}
  par(mfrow=c(2,3))
  plot(land_size,main="Population sizes")
  #plot(raster(migrationM))
  plot(raster(transitionmatrice),main="Transition matrix")
  plot(raster(a_value_simul),main="Simul genet differ")
  plot(raster(a_value_theory_stepping_stone_model),main="Expected stepping stone")
  plot(raster(a_value_theory_island_model),main="Expected island")
  plot(raster(a_value_theory_graph_model),main="Expected graph model")
  if (!genetDistis.null(file)){dev.off()}
  list(land_size,transitionmatrice,a_value_simul,a_value_theory_stepping_stone_model,a_value_theory_island_model,a_value_theory_graph_model,finalGenetData)
}



validation_with_coalescent <- function(rasterStack=rasterStack,
                                       pK = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),
                                                   nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),
                                                                               c("BIO1","BIO12"))), 
                                       pR = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),
                                                   nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),
                                                                               c("BIO1","BIO12"))),
                                       shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"),
                                       shapesR=c(BIO1="constant",BIO12="constant"), 
                                       shapeDisp="gaussian", 
                                       pDisp = .5, filen=NULL,
                                       nbLocus=60, initial_locus_value=200,stepvalue=2,
                                       mutation_model="stepwise",mutation_rate=1E-2,
                                       nind="Ne",
                                       stat="a_value",
                                       Option = "sample_1col_diploid")
{
  rasK=rasterStack;rasK[cellNumA(rasK)]= as.matrix(ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"]);rasK=rasK[[1]]
  geneticData = CreateGenetArray(rasK, nbLocus,initial_locus_value,Option=Option,nind=nind)
  #geneticData = CreateGenetArray(rasK, 4,200,Option="sample_1col_diploid",nind=8)
  coal_genet <- simul_coal_genet(geneticData,rasterStack,pK=pK,pR=pR,shapesK=shapesK,shapesR=shapesR,
                                 shapeDisp=shapeDisp,pDisp=pDisp,
                                 mutation_rate=mutation_rate,
                                 initial_locus_value=initial_locus_value,
                                 mutation_model=mutation_model,
                                 stepvalue=stepvalue,locusnames=NA)
  a_valueSimul <- a_value_matrix(coal_genet$tip_genotype,rasK) 
  FstSimul <- FstatsRas(coal_genet$tip_genotype,by="cell",all=TRUE,cells=cellNumA(rasterStack))$Fst
  genetDistSimul <- genetDist(tip_genotype = coal_genet$tip_genotype,method = "deltaMu",stepvalue = 2,byCell = TRUE)
      #a_value_theory_stepping_stone_model = distanceMatrix(rasterStack)/(4*0.05)
  if (shapeDisp %in% c("contiguous","gaussian")){
    stepping_migration_rate <- switch(shapeDisp,
                                      contiguous=pDisp,
                                      gaussian=2*dnorm(res(rasK)[1],0,pDisp))
    genetDistStepGM <- genetDistStepping(migrationMatrix(rasK,shapeDisp,pDisp),popSize = rasK,mutation_rate = mutation_rate)
    genetDistIslandGM <- matrix(1/(1+4*stepping_migration_rate*mean(valuesA(rasK)))/(1-1/(1+4*stepping_migration_rate*mean(valuesA(rasK)))),nrow=ncellA(rasK),ncol=ncellA(rasK))
    diag(genetDistIslandGM) <- 1/(1+4*(1-stepping_migration_rate*mean(valuesA(rasK))))
    if (any(dim(rasterStack) == ncell(rasK))) {#one dimensional landscape
      linearizedFst_theory_stepping_stone_model = distanceMatrix(rasK)/(res(rasK)[1])
      genetDistStepM <- 1/(stepping_migration_rate/2*distanceMatrix(rasK)/res(rasK))^2
    } else {#bidimensional landscape
      linearizedFst_theory_stepping_stone_model = log(distanceMatrix(rasterStack))/(res(rasterStack)[1])
      genetDistStepM <- log(1/(stepping_migration_rate/4*distanceMatrix(rasK)/res(rasK))^2)
    }
  } else {
    linearizedFst_theory_stepping_stone_model=matrix(NA,nrow=dim(transition)[1],ncol=dim(transition)[2])
    genetDistStepM=NA 
    genetDistStepMG=NA
  }
  linearizedFst_theory_island_model <- matrix(pDisp[1],nrow=ncellA(rasK),ncol=ncellA(rasK))-pDisp[1]*diag(ncellA(rasK))
  linearizedFst_theory_undigraph_model <- expect_linearizedFst_undigraph(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA)
  genetDist_Undigraph <- genetDistUndigraph(coal_genet$transition$backw,popSize=rasK,mutation_rate = mutation_rate)
  genetDist_Digraph <- genetDistDigraph(coal_genet$transition$backw,popSize=rasK,mutation_rate = mutation_rate)
  linearizedFst_theory_digraph_model <- expect_linearizedFst_digraph(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA)
  coalescence_times <- coalescence_prob_time_distribution_matrix(transition = coal_genet$transitionAndK$backw,max_time_interval = 4,rasK = rasK,threshold = 1E-5)
  coalescence_distance <- coalescence_times[["exp_times"]]*mutation_rate
  hitting_time <- hitting_time_digraph(coal_genet$transitionAndK$backw)
  sampleSize<- rasK
  cellsSampled <- as.numeric(levels(as.factor(geneticData$Cell_numbers)))
  for (i in cellNumA(rasK)) {values(sampleSize)[i] <- sum(coal_genet$tip_genotype$Cell_numbers==i)}
  genetDist_absorbingMethod <- genetDistAbsorbingMethod(coal_genet$transitionAndK$backw,valuesA(rasK),mutation_rate) 
  rasK<-as.matrix(rasK)
  
  result = list(coalescents_simulated=coal_genet$coalescent,
                tip_genotypes=coal_genet$tip_genotype,
                rasK=rasK,
                popSize = as.matrix(rasK),
                sampleSize=as.matrix(sampleSize),
                cellsSampled=cellsSampled,
                transitionmatrice=coal_genet$transitionAndK$backw,
                  #                a_value_Ind=a_value_Ind,
                a_valueSimul=a_valueSimul,
                Fst_Simul=FstSimul,
                linearizedFstSimul=FstSimul/(1-FstSimul),
                geneticDistSimul = genetDistSimul,
                linearizedFst_theory_undigraph_model=linearizedFst_theory_undigraph_model,
                linearizedFst_theory_digraph_model=linearizedFst_theory_digraph_model,
                genetDistUndigraph=  genetDist_Undigraph,
                genetDistDigraph=  genetDist_Digraph,
                hitting_time=hitting_time,
                coalescence_times_distribution=coalescence_times[[1]],
                coalescence_times_mean=coalescence_times[[2]],
                coalModelMutDist=coalescence_distance,
                genetDist_absorbingMethod=genetDist_absorbingMethod,
                parameters = list(rasterStack=rasterStack,pK=pK,pR =pR,
                                  shapesK=shapesK,shapesR=shapesR, shapeDisp=shapeDisp, 
                                  pDisp = pDisp, filen=filen,mutation_rate=mutation_rate, nbLocus=nbLocus, 
                                  initial_locus_value=initial_locus_value,stepvalue=stepvalue,
                                  mutation_model=mutation_model,mutation_rate=mutation_rate,
                                  nind=nind,stat=stat,Option=Option)
  )
  if (shapeDisp %in% c("contiguous","gaussian")){result[["linearizedFst_theory_stepping_stone_model"]]=linearizedFst_theory_stepping_stone_model
                                                 result[["genetDistStepping"]]=genetDistStepM
                                                 result[["genetDistGraphStepping"]]=genetDistStepGM
                                                 result[["genetDistIsland"]]=genetDistIslandGM
  }
result
}

#
#
# Analysis of landscape genetic data
#
#

samplePrior <- function(prior,method="random")
{
  p <- NA
  if (method=="random") {
    for (i in names(prior)){
      p[i] = switch(prior[[i]]$distribution,
                    uniform = runif(1,min=prior[[i]]$p["min"],max=prior[[i]]$p["max"]),
                    loguniform = exp(runif(1,min=log(prior[[i]]$p["min"]),max=log(prior[[i]]$p["max"]))),
                    fixed = prior[[i]]$p)
    }
  }
  if (method=="mean") {
    for (i in names(prior)){
      p[i] = switch(prior[[i]]$distribution,
                    uniform = mean(c(prior[[i]]$p["min"],prior[[i]]$p["max"])),
                    loguniform = exp(mean(c(min=log(prior[[i]]$p["min"]),max=log(prior[[i]]$p["max"])))),
                    fixed = prior[[i]]$p)
    }
  }
  p[-1]
}

landGenetAnalysis <- function(genetData,environmentalData,priors)
{
  # prior list is strictured as follows:
  # 4 sublists representing the parameter : $K, $R, $mutation_rate, $pDisp
  # each parameter has 3 components: $distribution is the name of the prior distribution 
  #                                  $model is the fuction or the model using the parameters
  #                                  $p is the prior hyper parameters
  # prior=list()
  # prior$K$distribution = "uniform"
  # prior$K$model = "proportional"
  # prior$K$p = c(min=0.001,max=0.5)
  # prior$R$distribution = "fixed"
  # prior$R$model  = "constant"
  # prior$R$p = 20
  # prior$mutation_rate$model = "stepwise"
  # prior$mutation_rate$distribution = "uniform"
  # prior$mutation_rate$p = array(c(0,1),dim=c(1,2),dimnames=list("mu",c("min","max")))
  # prior$pDisp$model="contiguous"
  # mrior$pDisp$distribution="uniform"
  # prior$pDisp$p=c(min=0.001,max=0.5)
  # prior = c(prior$pK,prior$ShapeK,prior$pR,prior$ShapeR,prior$mutation_model,prior$prior_mutation_rate,prior$ShapeDisp,prior$pDisp
  #  
  # NLM (error minimisationh, newto type algorithm)
  p=samplePrior(prior,method="random")
  gDist <- genetDist(genetData,method="deltaMu",stepvalue)
  genetDistUndigraphForNLM(p)
  parameters
}


#
# Function to check that we can jump generations in the coalescent time probability calculation
# checkTimeIntervall : produces all the expected coalescence times depending on the number of 
#                      generation jumped
# plot_coal_time_depending_on_time_interval : plots the results, 
#                                             and produces a table of correlations

checkTimeInterval <- function(transition = result$transitionmatrice,rasK = rasK,threshold = 1E-5)
{
  meanCoalTimes <- list()
  for (time_interval in c("1","4","10","15","20","40","100")){
    meanCoalTimes[[time_interval]] <- coalescence_prob_time_distribution_matrix(transition = transition,max_time_interval = as.numeric(time_interval),rasK = rasK,threshold = 1E-5)$exp_times
  }
  meanCoalTimes
}

#
# Function to check that we can jump generations in the coalescent time probability calculation
#

plot_coal_time_depending_on_time_interval <- function(meanCoalTimes,filen=NA,corr=TRUE,land=TRUE)
{
  if (land){
    if (!is.na(filen)) {
      if (grepl(".jpg",filen)) jpeg(paste("land_",filen,sep="")) else if (grepl(".pdf",filen)) pdf(paste("land_",filen,sep=""))
    }
    par(mfrow=c(3,3))
    for (i in 1:length(meanCoalTimes)){
      plot(raster(meanCoalTimes[[i]]),main=paste("jump",names(meanCoalTimes)[i]))
    }
    if (!is.na(filen)) dev.off()
  }
  if (corr){
    df<-as.data.frame(matrix(unlist(meanCoalTimes),nrow=prod(dim(meanCoalTimes[[1]])),ncol=length(meanCoalTimes),dimnames=list(1:prod(dim(meanCoalTimes[[1]])),names(meanCoalTimes))))
    if (!is.na(filen)) {
      if (grepl(".jpg",filen)|grepl(".jpeg",filen)) jpeg(paste("corr_",filen,sep="")) 
      if (grepl(".pdf",filen)) pdf(paste("corr_",filen,sep=""))
    }
    plot(df,cex=.1)    
    if (!is.na(filen)) dev.off()
    if (!is.na(filen)) {
      if (grepl(".jpg",filen)) jpeg(paste("corr_",filen,sep="")) 
      if (grepl(".pdf",filen)) pdf(paste("corr_",filen,sep=""))
    }
    if (!is.na(filen)) write.table(cor(df),file=paste("corr_",strsplit(filen,"\\.")[[1]][1],".csv",sep=""))    
    if (!is.na(filen)) dev.off()
  }
}
  
  
#
# Function to plot backward validation results
#


plot_validation <- function(validationdata,what=c("pop_size","sample_size","transition","geneticDistSimul",
                                                  "hitting_time","genetDistUndigraph",
                                                  "genetiDistDigraph","graphCoalTimes","graphCoalModelMutDist",
                                                  "Fsimul","FlinearSimul",
                                                  "linFstepping","genetDistStepping",
                                                  "genetDistIsland","FlinearDigraph",
                                                  "FlinearUndigraph","genetDistGraphStepping",
                                                  "genetDistAbsorbingMC"),
                            whichCell="all",
                            filen=NA)
{
  refWhat <- c("pop_size","sample_size","transition","geneticDistSimul",
               "hitting_time","genetDistUndigraph",
               "genetiDistDigraph","graphCoalTimes","graphCoalModelMutDist",
               "Fsimul","FlinearSimul",
               "linFstepping","genetDistStepping","genetDistIsland","FlinearDigraph",
               "FlinearUndigraph","genetDistGraphStepping","genetDistAbsorbingMC")
  if (is.numeric(what)) what <- refWhat[what]
  par(mar=c(3,3,3,3))
  cells<-as.numeric(levels(as.factor(validationdata$tip_genotypes$Cell_numbers)))
  n1=ceiling(sqrt(length(what)))
  n2=ceiling(length(what)/n1)
  if (!is.na(filen)) {ncols=n2;nlines=n1
                      if (grepl(".jp",filen)) jpeg(filen, width = 480*3, height = 480*3,res=300)
                      if (grepl(".pdf",filen)) pdf(filen)
                      } else {ncols=n1;nlines=n2}
  par(mfrow=c(nlines,ncols))
  if (is.numeric(whichCell)) {cells=whichCell} else if (is.character(whichCell)) {
    if (whichCell=="all") {cells=cellNumA(validationdata$rasK)}
    if (whichCell=="cells_sampled") {cells<-as.numeric(levels(as.factor(validationdata$tip_genotypes$Cell_numbers)))}
    if (grepl("KmoreThan",whichCell)) {cells=which(t(validationdata$rasK)>strsplit(whichCell,"Than")[[1]][2])}
    if (whichCell=="KmoreThan1") {cells=which(t(validationdata$rasK)>1)}
  } else {cells=cellNumA(validationdata$rasK)}
  title =c(pop_size="Populations sizes",
           sample_size="Samples sizes",
           transition = "transition",
           Fsimul ="F simul",
           FlinearSimul="F/(1-F) Simul",
           geneticDistSimul="Simulations",
           linFstepping="F/(1-F) Step",
           genetDistStepping="g dist step",
           genetDistIsland="Island",
           hitting_time = "Hitting time",
           FlinearDigraph = "F/(1-F) Digraph",
           FlinearUndigraph="F/(1-F) Undigraph",
           genetDistUndigraph="Undigraph",
           genetiDistDigraph="Digraph",
           graphCoalTimes="coalescence time",
           graphCoalModelMutDist="Markov chain",
           genetDistGraphStepping="Stepping Stone",
           genetDistAbsorbingMC="Absorbant MC")
  matrixNames <- c(pop_size="popSize",
                sample_size="sampleSize",
                transition = "transitionmatrice",
                Fsimul ="Fst_Simul",
                FlinearSimul="linearizedFstSimul",
                geneticDistSimul="geneticDistSimul",
                linFstepping="linearizedFst_theory_stepping_stone_model",
                genetDistStepping="genetDistStepping",
                hitting_time = "hitting_time",
                FlinearDigraph = "linearizedFst_theory_digraph_model",
                FlinearUndigraph="linearizedFst_theory_undigraph_model",
                genetDistUndigraph="genetDistUndigraph",
                genetiDistDigraph="genetDistDigraph",
                graphCoalTimes="coalescence_times_mean",
                graphCoalModelMutDist="coalModelMutDist",
                genetDistGraphStepping="genetDistGraphStepping",
                genetDistIsland="genetDistIsland",
                genetDistAbsorbingMC="genetDist_absorbingMethod")
  Cells1 <-  list(pop_size=1:dim(validationdata$rasK)[1],
              sample_size=1:dim(validationdata$rasK)[1],
              transition = 1:dim(validationdata$transitionmatrice)[1],
              Fsimul =cells,
              FlinearSimul=cells,
              geneticDistSimul=1:dim(validationdata$geneticDistSimul)[1],
              linFstepping=cells,
              genetDistStepping=cells,
              hitting_time = cells,
              FlinearDigraph = cells,
              FlinearUndigraph=cells,
              genetDistUndigraph=cells,
              genetiDistDigraph=cells,
              graphCoalTimes=cells,
              graphCoalModelMutDist=cells,
              genetDistGraphStepping=cells,
              genetDistAbsorbingMC=cells,
              genetDistIsland=cells)
  Cells2=Cells1
  Cells2[["pop_size"]] <- 1:dim(validationdata$rasK)[2]
  Cells2[["sample_size"]] <- 1:dim(validationdata$rasK)[2]
  for (mat in what){
    plot(raster(validationdata[[matrixNames[mat]]][Cells1[[mat]],Cells2[[mat]]]),main=title[mat])
  }
  if(!is.na(filen)) dev.off()
}



#test nlm
ssr <- function(p)
{
  popSize = K_Function(rasterStack = donneesEnvironmentObs, p[1], p[2:(nblayers+1)])
  matrice_migration = migrationMatrix(donneesEnvironmentObs,shapeDisp, p[(nblayers+2):(nbpDisp+nblayers+1)])
  matrice_transition = transitionMatrixBackward(popSize, matrice_migration)
  matrice_laplace = laplaceMatrix(matrice_transition)
  commute_time = resistDist(matrice_laplace)
  linearizedFst_att = linearizedFstUndigraph(commute_time, p[(nbpDisp+nblayers+2)])
  linearizedFst_att = linearizedFst_att[finalGenetData$Cell_Number,finalGenetData$Cell_Number]
  return(mean((linearizedFst_obs - linearizedFst_att)^2))
}

##
## Likelihood function
##
# transition =result$transitionmatrice; gDat=resutl$tip_genotype

formatGeneticData <- function(gDat,rasK)
{
  loci<- colnames(gDat)[grep("ocus",colnames(gDat))]
  ProbCellPairs<- array(0,dim=c(dim(gDat)[1]/2,length(loci),2,ncellA(rasK)),dimnames=list(paste("pair",1:(dim(gDat)[1]/2),sep=""),loci,c("gene1","gene2"),cellNumA(rasK)))
  GenetDistPairs <- array(NA,dim=c(dim(gDat)[1]/2,length(loci)),dimnames=list(paste("pair",1:(dim(gDat)[1]/2),sep=""),loci))
  for (locus in loci){
    gene1 <- sample(dim(gDat)[1],dim(Pairs)[1])
    gene2 <- sample((1:dim(gDat)[1])[-gene1])
    ProbCellPairs[,locus,1,gDat[gene1,"Cell_numbers"]] <- 1
    ProbCellPairs[,locus,2,gDat[gene2,"Cell_numbers"]] <- 1
    GenetDistPairs[,locus] <- abs(gDat[gene1,locus]-gDat[gene2,locus])
  }
list(ProbCellPairs=ProbPairs,GenetDistPairs=GenetDistPairs)
}

deme_coocurence_probability <- function(pGenes,transition,time)
{ # pGenes is a matrix of initial occurence  probabilities for a set of genes
  # transition is a transiiton matrix per unit of time
  # time is time
  # t(M^t P0[g1])%*%(M^t P0[g2])
# P0G1
  p=pGene1%*%t(pGene2)
}

likelihoodCoalescent <- function(coalescent)
{
  
}

probgenet <- function(transition,gDat)
{
  if (all(unlist(strsplit(colnames(gDat)[grep("ocus",colnames(gDat))],"\\."))
          [1:length(colnames(gDat)[grep("ocus",colnames(gDat))])*2]==
            rep(c("1","2"),length(colnames(gDat)[grep("ocus",colnames(gDat))])))){# gDat is Two cols genotype format
    gDat <- TwoCols2OneCol(gDat)
  } 
  for (allele in grep("ocus",colnames(gData))){
    
  }
}


##
##
## Functions unused at the moment
##


S <- function(rasterStack)
{# S from Slatkin 1993 Evolution Eq. 8b
  S=matrix(0,nrow=ncellA(rasterStack),ncol=ncellA(rasterStack))
  d=dim(rasterStack)[1:2]
  # we get the lines and column number of each cell of sthe rasterStack
  ij <- t(t(xyFromCellA(rasterStack))/(res(rasterStack))+.5)
  for (donor in 1:ncellA(rasterStack)){
    for (receiver in 1:ncellA(rasterStack)){
      i <- abs(ij[donor,1]-ij[receiver,1])
      j <- abs(ij[donor,2]-ij[receiver,2])
      for (k in 0:(d[1]-1)){
        for (l in 0:(d[2]-1)){
          S[donor,receiver]=S[donor,receiver]+(2-as.numeric(k==0))*(2-as.numeric(l==0))*(1-cos(pi*i*k/d[1])*cos(pi*j*l/d[2]))/(1-.5*(cos(pi*k/d[1])+cos(pi*l/d[2])))
          # modified from slatkin's torus: considering the landscape is represented by a torus of double number of demes in each direction
        }
      }
    }
  }
  S
}
