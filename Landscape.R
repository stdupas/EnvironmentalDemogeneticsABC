library(raster)
##########SET CLASS#########################################

setClass("Landscape",
         contains = "RasterStack",
         slots = c(period="Date",vars="character"),
         validity = validLandscape
         )
# j'ai retiré : Landscape<- car tu définis la fonction après
# j'ai modifié la longueur d la période ne peut être != 2
# jai fait d'autres checkin

validLandscape = function(object){
           if (length(object@period)!=2) stop("the period is  not valid because it contains more or less than two dates")
           if (object@period[2]<object@period[1]) stop("the period is not valid because the starting later than the ending")
           if (length(dim(object))!=3) stop("error when creating landscape : the array does not have 3 dimensions")
           if(length(object@vars)!=dim(object)[3]) stop("error when creating landscape : the number of layers in the array differs from the number of variables in var")
           if(any(dimnames(object)[[3]]!=object@vars)) stop("error when creating landscape : names of third dimension of landscape array do not correspond to slot vars names")
         TRUE
         }


Landscape<-function(rasterstack=a,period=p, vars=l){
  if (length(period)==1)  period<-c(period,period)
  names(rasterstack) <- vars
  new("Landscape",rasterstack,period=period,vars=vars)
}

## tu as pas besion de méthode is.landscape (cf plus loin)
## il suffit de mettre class(object) pour savoir si c'est un landscape
# j'ai ajouté les noms de var pour les variables de l'array
# attention length(period)/2 peut être impair
# j'ai retiré :
#       if(period[1]!=period[2*i-1]&&period[2]!=period[2*i])stop("The period is not the same for each layer")
# en fait on ne va proposer qu'une seule  période pour tout l'objet
# une période est un vecteur de dates de longueur 2 pour laquelle la seconde  date est supérieure à la première

LandscapeHistory<-setClass("LandscapeHistory",
                           contains = "list",
                           validity = function(object){
                             if(any(unlist(lapply(1:length(object),function(x) class(object[[x]])!="Landscape")))) stop("An element of the list is not a Landscape")
                             if (any(unlist(lapply(1:length(object),function(x) any(object[[x]]@vars!=object[[1]]@vars))))) stop("error in lanscape list, vars differ between periods")
                           }
                           )


# j'ai modifié is.landscape(object) car l'objet n'est pas un landscape, mais une liste!
# j'ai mis à la place un lapply 
# j'ai viré LandscapeHistory<-function(a){ car la fonciton était déja définie
## attention l'objet est une  liste, ce ne sera pas un landscape
# is. landscape va pas marcher
# il faut toujour etre très explicite  dans les texte d'erreur. Fais plus attention à la forme.

nbpar <- function(x) {unlist(lapply(x,function(x) switch(x,
                            constant=1,
                            proportional=1,
                            enveloppe=2,
                            envelin=2,
                            quadratic=4)
                            ))
                      }


setClass("NicheModel",
         slots = c(variables="character",parameterList="list",reactNorms="character",form="character"),
         validity=validityNicheModel
)

setClass("form",
		slots=c(numberVar="integer",additions="integer",multiplications="integer",openParenthesis="integer",closeParenthesis="integer"),
      validity=function(object){
        if (length(object@openPatrenthesis)!=length(object@closeParenthesis)) stop("error in 'form': different number of openParenthesis and closeparenthesis")
        if (numberVar<=length(additions)+length(multiplications)) stop("error in 'form': too much operations")

        if (order(append(addition,multiplications))) stop("error in 'form': too much operations")
  }
)


NicheModel<-function(variables=var,parameterList=par,reactNorms=rea,period=pe,form=formul){
  names(parameterList)=variables
  names(reactNorms)=variables
  new("NicheModel",variables=variables,parameterList=parameterList,reactNorms=reactNorms,form=form)
}

validityNicheModel = function(object){
              if(FALSE%in%lapply(object@parameterList,is.numeric))stop("error in NicheModel parameter list : Parameter just accept numeric!")
              if(length(object@variables)!=length(object@reactNorms))stop("error in NicheModel : number of variables and number of reaction norms do not correspond")
              notMatching <- (nbpar(object@reactNorms) != unlist(lapply(1:length(object@parameterList),function(x) length(object@parameterList[[x]])))) 
              if (any(notMatching)) stop(paste("error in NicheModel : number of paremeters and reactionNorm do not match for variable ",which(notMatching),". ",sep=""))
#              if grep("(",object@form)
TRUE 
}


#######SET METHODS##########################################


setGeneric(
  name = "runNicheModel",
  def=function(object,model){return(standardGeneric("runNicheModel"))}
)



setMethod("runNicheModel",
          signature=c("Landscape","NicheModel"),
          definition = function(object,model){                  #X=object, p=,shape=
            Y=lapply(model@variables,function(x){
                   									switch(model@reactNorms[[x]],
                   									       constant={object[[x]]=model@parameterList[[x]]},
                   									       proportional = {object[[x]]=as.matrix(object[[x]])*model@parameterList[[x]]},
                   									       enveloppe = {object[[x]]=envelope(as.matrix(object[[x]]),model@parameterList[[x]])},
                   									       envelin={object[[x]]=envelinear(as.matrix(object[[x]]),model@parameterList[[x]])},
                   									       conQuadratic={object[[x]]=conQuadratic(as.matrix(object[[x]]),model@parameterList[[x]])} 
                                                 #conquadraticskewed=conquadraticskewed(object[,,(model@variables==x)],p),
                                                 #conquadraticsq=conquadraticsq(object[,,(model@variables==x)],p),
                                                 #conquadraticskewedsq=conquadraticskewedsq(object[,,(model@variables==x)],p)
                   									)
                         }
                )
            Y=lapply(Y,prod)
          }
)


envelope <- function(X,p){
X>=p[1]&X<=p[2]
}

envelinear <- function(X, p) {
(X-p[1])/(p[2]-p[1])*envelope(X,p)
}
plot(1:100,envelinear(1:100,c(20,60)))

constant <- function(X,p){X[]<-p}

conQuadratic <- function(X,p)
{
-4*(X-p[2])*(X-p[1])/((p[2]-p[1])^2)*envelope(X,p)
}
plot(1:100,conQuadratic(1:100,c(20,60)))

conQuadratic <- function(X,p){
quadraticConcave(X,p)*envelinear(X,p)
}
plot(1:100,quadraticConcaveSkewed(1:100,c(20,60)))

#fonction quadratique bornée par p1 et p2 et dont le maximum est 1
#f(x)=-ax²+bx+c
#f(x)=0
#f(x)=-a(X-p1)*(X-p2)=-aX²+a(p1+p2)X+ap1p2 
#f(p1)=0;f(p2)=0
#df/dx =-2ax+a(p1+p2)
#df/dx((p1+p2)/2)=-2a
#1=a(p2-p1)/2*(p1-p2)/2
#a=-4/(p2-p1)^2

##########MANIPULATION CLASS################################

r <- raster(ncol=40, nrow=20)
r[] <- rnorm(n=ncell(r))
s <- stack(x=c(r, r*2, r+1,r*3))
per<-as.Date( c("2017-02-01","2017-02-01"))
pe1 <- as.Date(c("2007-01-01","2007-12-31"))
pe2 <- as.Date(c("2008-01-01","2009-12-31"))
vari<-c("l","t","p","h")
vari<-c("l","t","p")
para<-list(c(1,3),2,c(4,5.2))
para<-list(1,2,c(4,5.2),c(3,8))
rea<-c(l="constant",t="proportional",p="enveloppe",h="envelin")
rea<-c(l="constant",t="constant",p="constant")
formul=c(1,"*","(",2,"*",3,"*",4,")")



lscp1<-Landscape(rasterstack = s,period=as.Date("2017-02-01"),vars=vari)

lscp2<-Landscape(Array=array(12:1,dim=c(2,2,3)),period=as.Date(c("2017-02-02","2017-02-06")),vars=vari)


model<-NicheModel(variables=vari,parameterList=para,reactNorms=rea,form=formul)
landhistory <- LandscapeHistory(list(Landscape(LandscapeArray1,period=pe1,vars=vari),Landscape(LandscapeArray2,period=pe2,vars=vari)))
lista<-list(lscp1,lscp2)
lh1<-LandscapeHistory(lista)

a<-runNicheModel(lscp1,model)
a
