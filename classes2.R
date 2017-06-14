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



TransitionBackward<-setClass("TransitionBackward",
                             contains = "matrix",
                             prototype = prototype(matrix(nrow=100,ncol=100)),
                             validity = function(object){
                               if (nrow(object)==0)stop("The matrix is empty.")
                               if (nrow(object)!=ncol(object))stop("The matrix is not square")
                               if (all(rowSums(object)==1))TRUE else stop("The sum of probabilities in each row is not 1.")
                               }
                             )

 
#new(Class="TransitionBackward",matrix=mat)
TransitionBackward<- function(matrix){
  if (nrow(matrix)!=ncol(matrix))stop("The matrix is not square")
  if(class(rownames(matrix)[1])!="character"){
      lname <- c(1:nrow(matrix))
      rownames(matrix) <- lname
      colnames(matrix) <- lname
      }
  new(Class="TransitionBackward",matrix)
}

mat <- matrix(data=5:10, nrow=2, ncol=3, byrow=T)
mat2 <- matrix( nrow=0, ncol=0, byrow=T)
mat3 <- rbind(c(0.2,0.8),c(0.6,0.4))
mat4 <- rbind(c(0.2,0.7),c(0.3,0.4))

mic<-TransitionBackward(mat3)


laplaceMatrix(mic)
commute_time_undigraph(mic)
hitting_time_digraph(mic)
getMatrix(mic)
mic@matrix
