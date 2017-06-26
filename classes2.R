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


TransitionBackward<- function(matrix){
  if (nrow(matrix)!=ncol(matrix))stop("The matrix is not square")
  if(class(rownames(matrix)[1])!="character"){
      lname <- c(1:nrow(matrix))
      rownames(matrix) <- lname
      colnames(matrix) <- lname
      }
  new(Class="TransitionBackward",matrix)
}



CoalescentModel<-setClass("CoalescentModel",
                    contains = "Demographic",
                    prototype = prototype(TransitionBackward=TransitionBackward(mic)),
                    validity = function(object){
                      if (is.null(TransitionBackward))stop("The TransitionBackward is NULL")
                    }
                    )

Demographic<-setClass("Demographic",
                      slots = c(K="numeric", R="numeric"),
                      prototype = prototype(K=1,R=1),
                      validity = function(object){
                        if(K<0)stop("K is negative")
                        if(R<0)stop("R is negative")
                      }
                      )


Genetic<-setClass("Genetic",
                  prototype = prototype(geneticData=data.frame(1, 1:10, sample(LETTERS[1:3], 10, replace = TRUE))),
                  representation = representation(
                    geneticData = "data.frame"
                  ),
                  validity = function(object){
                    if (is.null(object@geneticData))stop("The dataBase is null.")
                  }
)


genetic <- function(df=data.frame(locus1.1=c(200,202),locus1.2=c(204,200)),ploidy=NA, ploidyByrow=NA){
  if (is.na(ploidy)) { 
    if (any(grep("\\.2", colnames(df)))) 
    {
      P <- (unlist(strsplit(grep("\\.",colnames(df),value=TRUE),"\\.")))
      ploidy <- max(suppressWarnings(as.integer(P))[!is.na(suppressWarnings(as.integer(P)))])
    }
    if (any(grep("\\.2", rownames(df)))) 
    {
      P <- (unlist(strsplit(grep("\\.",rownames(df),value=TRUE),"\\.")))
      ploidy <- max(suppressWarnings(as.integer(P))[!is.na(suppressWarnings(as.integer(P)))])
    }
  }
  if (is.na(ploidyByrow)) ploidyByrow = !(any(grep(paste("\\.",ploidy,sep=""), colnames(df))))
  new("genetic",df,ploidy=ploidy,ploidyByrow=ploidyByrow)
}


spatialGenetic <- setClass("spatialGenetic",
                           slots = c(x="numeric", y="numeric",Deme_numbers="integer"),
                           contains = "genetic",
                           prototype = prototype(genetic(),x=c(1,2),y=c(1,1),Deme_numbers=as.integer(c(1,2))),
                           validity = function(object){
                             if (length(object@x)!=length(object@y)) stop("slots x and y do not have the same length")
                             if (length(object@x)!=nrow(object)) stop("slots x and genetic do not have the same number of individuals")
                           }
)


spatialGenetic <- function(df=NA,x=NA,y=NA,Deme_numbers=NA)
{
  if(is.na(df)) df=cbind(data.frame(genetic()),x=c(1,2),y=c(1,1),Cell_numbers=c(1,2))
  if (!(all(c("x","y","Deme_numbers")%in%colnames(df)))){
    if (is.na(Deme_numbers)){
      if ("Deme_numbers"%in%colnames(df)) {
        Deme_numbers=df$Deme_numbers
        df <- df[,-which(colnames(df)=="Deme_numbers")]
      }
    }
    if (is.na(x)){
      if ("x"%in%colnames(df)) {
        x=df$x
        df <- df[,-which(colnames(df)=="x")]
      }
    }
    if (is.na(y)){
      if ("x"%in%colnames(df)) {
        y=df$y
        df <- df[,-which(colnames(df)=="y")]
      }
    }
  }
  new("spatialGenetic",genetic(df[,-which(colnames(df)%in%c("x","y","Deme_numbers"))]),x=df$x,y=df$y,Deme_numbers=df$Deme_numbers)
}




d1<-data.frame()






mat <- matrix(data=5:10, nrow=2, ncol=3, byrow=T)
mat2 <- matrix( nrow=0, ncol=0, byrow=T)
mat3 <- rbind(c(0.2,0.8),c(0.6,0.4))
mat4 <- rbind(c(0.2,0.7),c(0.3,0.4))

mic<-TransitionBackward(mat3)
mic["matrix"]
mic

laplaceMatrix(mic)
commute_time_undigraph(mic)
hitting_time_digraph(mic)
genetDistUndigraph(mic,100,0.1)
genetDistDigraph(mic,100,0.1,method = "Goldstein95")

lapply()
matrix