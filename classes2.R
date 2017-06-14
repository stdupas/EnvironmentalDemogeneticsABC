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



TransitionBackward<-setClass("TansitionBackward",
                             contains = "matrix",
                             prototype = prototype(matrix(nrow=100,ncol=100)),
                             validity = function(object){
                               if (nrow(object@matrix)==0)stop("The matrix is empty.")
                               if (all(rowSums(object@matrix)==1))TRUE else stop("The sum of probabilities in each row is not 1.")
                               }
                             )

 
#new(Class="TransitionBackward",matrix=mat)
#TransitionBackward<- function(matrix=matrix){
#  new(Class="TransitionBackward",matrix=matrix)
#}

mat <- matrix(data=5:10, nrow=2, ncol=3, byrow=T)
mat2 <- matrix( nrow=0, ncol=0, byrow=T)
mat3<-rbind(c(0.2,0.5,0.3),c(0.6,0.1,0.3))

mic<-TransitionBackward(matrix = mat3)

laplaceMatrix(mic)
getMatrix(mic)
mic@matrix
