library(raster)
##########SET CLASS#########################################

setClass("Landscape",
         contains = "RasterStack",
         slots = c(period="Date",vars="character",distanceMatrix="matrix"),
         validity = validLandscape
         )

validLandscape = function(object){
           if (length(object@period)!=2) stop("the period is  not valid because it contains more or less than two dates")
           if (object@period[2]<object@period[1]) stop("the period is not valid because the starting later than the ending")
           if (length(dim(object))!=3) stop("error when creating landscape : the array does not have 3 dimensions")
           if(length(object@vars)!=dim(object)[3]) stop("error when creating landscape : the number of layers in the array differs from the number of variables in var")
           if(any(names(object)!=object@vars)) stop("error when creating landscape : names of third dimension of landscape array do not correspond to slot vars names")
           TRUE
         }


Landscape<-function(rasterstack=a,period=p, vars=l){
  if (length(period)==1)  period<-c(period,period)
  names(rasterstack) <- vars
  b<-xyFromCellA(rasterstack)
  mat=sapply(1:nrow(b),function(l1){
    sapply(1:nrow(b),function(l2){
      sqrt((b[l1,1]-b[l2,1])^2+(b[l1,2]-b[l2,2])^2)
    })})
  new("Landscape",rasterstack,period=period,vars=vars,distanceMatrix=mat)
}


LandscapeHistory<-setClass("LandscapeHistory",
                           contains = "list",
                           validity = function(object){
                             if(any(unlist(lapply(1:length(object),function(x) class(object[[x]])!="Landscape")))) stop("An element of the list is not a Landscape")
                             if (any(unlist(lapply(1:length(object),function(x) any(object[[x]]@vars!=object[[1]]@vars))))) stop("error in lanscape list, vars differ between periods")
                           }
                  )


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


setClass("MigrationModel",
         slots = c(shapeDisp="ANY",pDisp="ANY"),
         validity = validityMigrationModel
)

validityMigrationModel=function(object){
  if(!is.character(object@shapeDisp))stop("error in  MigrationModel shapeDisp : Parameter just accept character!")
  if(FALSE%in%lapply(object@pDisp,is.numeric))stop("error in MigrationModel pDisp : pDisp just accept numeric!")
  TRUE
}

MigrationModel<-function(shape=s,param=p){
  new("MigrationModel",shapeDisp=shape,pDisp=param)
}

setClass("EnvDinModel",
         slots = c(K="NicheModel",R="NicheModel",migration="MigrationModel"),
         validity=function(object){
           if(class(object@K)!="NicheModel")stop("Error in envDinModel K : K just accept NicheModel!")
           if(class(object@R)!="NicheModel")stop("Error in envDinModel R : R just accept NicheModel!")
           if(class(object@migration)!="MigrationModel")stop("Error in envDinModel migration : migration just accept MigrationModel!")
         }
         )

EnvDinModel<-function(K=k,R=r,migration=m){
  new("EnvDinModel",K=K,R=R,migration=migration)
}


#######SET METHODS##########################################

setMethod(
  f ="[",
  signature = c(x="EnvDinModel" ,i="character",j="missing"),
  definition = function (x ,i ,j , drop ){
    switch ( EXPR =i,
             "K" ={return(x@K)} ,
             "R" ={return(x@R)} ,
             "migration" ={return(x@migration)} ,
             stop("This slots doesn't exist!")
    )
  }
)

setMethod(
  f ="[",
  signature = c(x="MigrationModel" ,i="character",j="missing"),
  definition = function (x ,i ,j , drop ){
    switch ( EXPR =i,
             "pDisp" ={return(x@pDisp)} ,
             "shapeDisp" ={return(x@shapeDisp)} ,
             stop("This slots doesn't exist!")
    )
  }
)

setMethod(
  f ="[",
  signature = c(x="Landscape" ,i="character",j="missing"),
  definition = function (x ,i ,j , drop ){
    switch ( EXPR =i,
             "period" ={return(x@period)} ,
             "distanceMatrix" ={return(x@distanceMatrix)} ,
             "vars" ={return(x@vars)} ,
             stop("This slots doesn't exist!")
    )
  }
)


setGeneric(
  name = "runNicheModel",
  def=function(object,model){return(standardGeneric("runNicheModel"))}
)



setMethod("runNicheModel",
          signature=c("Landscape","NicheModel"),
          definition = function(object,model){                  #X=object, p=,shape=
            Y=lapply(model@variables,function(x){
                   									switch(model@reactNorms[[x]],
                   									       constant={object[[x]]<-setValues(object[[x]],rep(model@parameterList[[x]],ncell(object[[x]])))},
                   									       #proportional = {values(object[[x]])=object[[x]]*model@parameterList[[x]]},
                   									       enveloppe = {object[[x]]=envelope(object[[x]],model@parameterList[[x]])},
                   									       envelin={object[[x]]=envelinear(object[[x]],model@parameterList[[x]])},
                   									       conQuadratic={object[[x]]=conQuadratic(object[[x]],model@parameterList[[x]])} 
                                                 #conquadraticskewed=conquadraticskewed(object[,,(model@variables==x)],p),
                                                 #conquadraticsq=conquadraticsq(object[,,(model@variables==x)],p),
                                                 #conquadraticskewedsq=conquadraticskewedsq(object[,,(model@variables==x)],p)
                   									)
                         }
                )
            Y=prod(stack(Y))
          }
)



envelope <- function(X,p){
  if(length(p)!=2)stop("The parameter is  not valid because it contains more or less than two values")
  else X>=p[1]&X<=p[2]
}

envelinear <- function(X, p) {
  if(length(p)!=2)stop("The parameter is  not valid because it contains more or less than two values")
  else (X-p[1])/(p[2]-p[1])*envelope(X,p)
}
#plot(1:100,envelinear(1:100,c(20,60)))

constant <- function(X,p){X[]<-p}

conQuadratic <- function(X,p)
{
  if(length(p)!=2)stop("The parameter is  not valid because it contains more or less than two values")
  else -4*(X-p[2])*(X-p[1])/((p[2]-p[1])^2)*envelope(X,p)
}
#plot(1:100,conQuadratic(1:100,c(20,60)))

conQuadraticsKed <- function(X,p){
quadraticConcave(X,p)*envelinear(X,p)
}
#plot(1:100,quadraticConcaveSkewed(1:100,c(20,60)))

#fonction quadratique bornée par p1 et p2 et dont le maximum est 1
#f(x)=-ax²+bx+c
#f(x)=0
#f(x)=-a(X-p1)*(X-p2)=-aX²+a(p1+p2)X+ap1p2 
#f(p1)=0;f(p2)=0
#df/dx =-2ax+a(p1+p2)
#df/dx((p1+p2)/2)=-2a
#1=a(p2-p1)/2*(p1-p2)/2
#a=-4/(p2-p1)^2
#########CREER TRANSITION MATRIX#############################

setGeneric(
  name = "createTransitionMatrix",
  def=function(object,model){return(standardGeneric("createTransitionMatrix"))}
)

setMethod(f="createTransitionMatrix",
          signature=c("Landscape","EnvDinModel"),
          definition=function(object,model){
            lpar<-runEnvDinModel(object,model)
            if ((length(lpar$R)==1)&(length(lpar$K)==1)){transition = lpar$R * lpar$K * t(lpar$migration)}
            if ((length(lpar$R)>1)&(length(lpar$K)==1)){transition = t(matrix(lpar$R,nrow=length(lpar$R),ncol=length(lpar$R))) * lpar$K * t(lpar$migration)}
            if ((length(lpar$R)==1)&(length(lpar$K)>1)){transition = lpar$R * t(matrix(lpar$K,nrow=length(lpar$K),ncol=length(lpar$K))) * t(lpar$migration)}
            if ((length(lpar$R)>1)&(length(lpar$K)==1)){transition = t(matrix(lpar$R,nrow=length(lpar$R),ncol=length(lpar$R))) * lpar$K * t(lpar$migration)}
            if ((length(lpar$R)>1)&(length(lpar$K)>1)) {transition = t(matrix(lpar$R,nrow=length(lpar$R),ncol=length(lpar$R))) * t(matrix(lpar$K,nrow=length(lpar$K),ncol=length(lpar$K))) * t(lpar$migration)}
            #TransitionBackward(transition)
          }
)



setGeneric(
  name = "runEnvDinModel",
  def=function(object,model){return(standardGeneric("runEnvDinModel"))}
)


setMethod(f="runEnvDinModel",
          signature=c("Landscape","EnvDinModel"),
          definition=function(object,model){
            R<-values(runNicheModel(object,model["R"]))
            K<-values(runNicheModel(object,model["K"]))
            migrationMat<-migrationMatrix(object,model["migration"])
            list(R=R,K=K,migration=migrationMat)
          }
)


setGeneric(
  name = "migrationMatrix",
  def=function(object,model){return(standardGeneric("migrationMatrix"))}
)



setMethod(
  f="migrationMatrix",
  signature=c("Landscape","MigrationModel"),
  definition=function(object,model)
  {
    Ndim = 1+all(ncell(object)!=dim(object)[1:2])
    migration = apply(object["distanceMatrix"], c(1,2), 
                      function(x)(switch(model["shapeDisp"],
                                         fat_tail1 = 1/(1+x^model["pDisp"][2]/model["pDisp"][1]),
                                         gaussian = (dnorm(x, mean = 0, sd = model["pDisp"][1], log = FALSE)),
                                         exponential = (dexp(x, rate = 1/model["pDisp"][1], log = FALSE)),
                                         contiguous = (x==0)*(1-model["pDisp"][1])+((x>0)-(x>1.4*res(object)[1]))*(model["pDisp"][1]/(2*Ndim)),
                                         contiguous8 = (x==0)*(1-object["pDisp"][1])+((x>0)-(x>2*res(object)[1]))*(model["pDisp"][1]/(4*Ndim)),
                                         island = (x==0)*(1-model["pDisp"][1])+(x>0)*(model["pDisp"][1]),
                                         fat_tail2 = x^model["pDisp"][2]*exp(-2*x/(model["pDisp"][1]^0.5)),
                                         contiguous_long_dist_mixt = model["pDisp"]["plongdist"]/ncellA(object)+(x==0)*(1-model["pDisp"]["pcontiguous"]-model["pDisp"]["plongdist"])+((x>0)-(x>1.4*res(object)[1]))*(model["pDisp"]["pcontiguous"]/2),
                                         gaussian_long_dist_mixt = model["pDisp"][2]/ncellA(object) + (dnorm(x, mean = 0, sd = model["pDisp"][1], log = FALSE))
                      )))
        return(migration/sapply(rowSums(migration),function(x)rep(x,ncol(migration))))
  }
)

setGeneric(
  name = "xyFromCellA",
  def=function(object){return(standardGeneric("xyFromCellA"))}
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterLayer",
  definition = function(object){
    df=xyFromCell(object,cellNumA(object))
    rownames(df) <- cellNumA(object)
  }
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterStack",
  definition = function(object){
    df =xyFromCell(object,cellNumA(object))
    rownames(df) <- cellNumA(object)
    df
  }
)

setGeneric(
  name = "nCellA",
  def=function(object){return(standardGeneric("nCellA"))}
)

setMethod(
  f = "nCellA",
  signature = "RasterLayer",
  definition = function(object){
    length(na.omit(values(object)))
  }
)

setMethod(
  f = "nCellA",
  signature = "RasterStack",
  definition = function(object){
    ncellA(object[[1]])
  }
)

setGeneric(
  name = "valuesA",
  def=function(object){return(standardGeneric("valuesA"))}
)

setMethod(
  f = "valuesA",
  signature = "RasterLayer",
  definition = function(object){
    #x=data.frame(variable=na.omit(values(object)))
    select <- !is.na(values(object))
    x=values(object)[select]
    names(x) <- which(select)
    #colnames(x)=names(object)
    x
  }
)

setMethod(
  f = "valuesA",
  signature = "RasterStack",
  definition = function(object){
    x=na.omit(values(object))
    colnames(x)=names(x)
    rownames(x) <- cellNumA(object)
    x
  }
)

setGeneric(
  name = "cellNumA",
  def=function(object){return(standardGeneric("cellNumA"))}
)

setMethod(
  f = "cellNumA",
  signature = "RasterLayer",
  definition = function(object){
    which(!is.na(values(object)))
  }
)

setMethod(
  f = "cellNumA",
  signature = "RasterStack",
  definition = function(object){
    cellNumA(object[[1]])
  }
)




##########MANIPULATION CLASS################################

r <- raster(ncol=2, nrow=2)
r[] <- c(c(1,2),c(3,4))#rnorm(n=ncell(r))
s <- stack(x=c(r, r*2, r+1,r*3))

class(r[])

vari<-c("l","t","p","h")
#vari<-c("l","t","p")
#para<-list(c(1,3),2,c(4,5.2))
para<-list(1,c(1,5),c(1,5.2),c(3,8))
rea<-c(l="constant",t="envelin",p="enveloppe",h="envelin")
#rea<-c(l="constant",t="constant",p="constant")
formul=c(1,"*","(",2,"*",3,"*",4,")")



lscp1<-Landscape(rasterstack = s,period=as.Date("2017-02-01"),vars=vari)
plot(lscp1)
for(i in 1:4){
  plot(lscp1[[i]])
}
lscp1
lscp2<-Landscape(Array=array(12:1,dim=c(2,2,3)),period=as.Date(c("2017-02-02","2017-02-06")),vars=vari)


model<-NicheModel(variables=vari,parameterList=para,reactNorms=rea,form=formul)

landhistory <- LandscapeHistory(list(Landscape(LandscapeArray1,period=pe1,vars=vari),Landscape(LandscapeArray2,period=pe2,vars=vari)))
lista<-list(lscp1,lscp2)
lh1<-LandscapeHistory(lista)

a<-runNicheModel(lscp1,model)
lscp1["distanceMatrix"]
class(lscp1[1])
values(a)
setValues()
class(values(a))
a
plot(a)
par(mfrow=c(2,3))
valuesA(a)

m<-MigrationModel(shape="gaussian",param = 1)
m["pDisp"]
edm1<-EnvDinModel(K=model,R=model,migration = m)
createTransitionMatrix(lscp1,edm1)
edm1
transi1<-createTransitionMatrix(lscp1,edm1)



###################################################################################"
r1<- raster(ncol=2, nrow=2)
r2<- raster(ncol=2, nrow=2)
r2[] <- rep(2,2:2)
r1[] <- rep(1,2:2)
s <- stack(x=c(r1,r2))
land1<-Landscape(rasterstack = s,period = as.Date(c("2017-02-02","2017-02-06")),vars = c("l","p"))
model1<-NicheModel(var=c("l","p"),reactNorms = c(l="constant",p="constant"),parameterList = list(0.5,1),form =c("(",1,"*",2,")") )
migra1<-MigrationModel(shape = "fat_tail1",param = c(1,1))
env1<-EnvDinModel(K=model1,R=model1,migration = migra1)
a=(runEnvDinModel(land1,env1))
a
b=createTransitionMatrix(land1,env1)
b





