library(raster)
######### PRECURSEUR ######################################################

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
######### SET CLASS ########################################################################################################

validLandscape = function(object){
  if (length(object@period)!=2) stop("the period is  not valid because it contains more or less than two dates")
  if (object@period[2]<object@period[1]) stop("the period is not valid because the starting later than the ending")
  TRUE
}

setClass("Landscape",
         contains = "RasterStack",
         slots = c(period="Date",vars="character",distanceMatrix="matrix"),
         validity = validLandscape
         )


Landscape<-function(rasterstack=rasterstack,period=dateVector, vars=charVector){
  if(class(rasterstack)!="RasterStack")stop("error in Landscape rasterstack : rasterstack just accept RasterStack!")
  if (length(rasterstack@layers)==0)stop("rasterstack values is null!")
  if(class(period)!="Date")stop("error in Landscape period : period just accept Date!")
  if(!is.character(vars))stop("error in Landscape vars : vars just accept character!")
  if (length(period)==1)  period<-c(period,period)
  if(length(names(rasterstack))!=length(vars)) stop("error when creating landscape : the number of layers in the rasterstack differs from the number of variables in var")
  names(rasterstack) <- vars
  b<-xyFromCellA(rasterstack)
  mat=sapply(1:nrow(b),function(l1){
    sapply(1:nrow(b),function(l2){
      sqrt((b[l1,1]-b[l2,1])^2+(b[l1,2]-b[l2,2])^2)
    })})
  new("Landscape",rasterstack,period=period,vars=vars,distanceMatrix=mat)
}


setClass("LandscapeHistory",
                           contains = "list",
                           validity = function(object){
                             if(!is.list(object))stop("error in lanscape list : landscape list just accept list.")
                             if(any(unlist(lapply(1:length(object),function(x) class(object[[x]])!="Landscape")))) stop("An element of the list is not a Landscape.")
                             if (any(unlist(lapply(1:length(object),function(x) any(object[[x]]@vars!=object[[1]]@vars))))) stop("error in lanscape list, vars differ between periods.")
                             if(length(object)>1){
                               if(any(unlist(lapply(1:length(object),function(x) lapply(x:length(object),function(y)if(x!=y)any(object[[y]]@period==object[[x]]@period))))))stop("error in lanscape period, at least two landscape have same period.")
                               if(any(unlist(lapply(1:(length(object)-1),function(x) (object[[x]]["period"][2]+1)!=object[[x+1]]["period"][1]))))stop("error in lanscape list, periods are not contigous")
                             }
                           }
                  )

LandscapeHistory<-function(Landscapelist=listOfLandscape){
  li<-unlist(lapply(1:length(Landscapelist),function(x) Landscapelist[[x]]["period"][1]))
  o<-order(c(li))
  lo<-unlist(lapply(1:length(o),function(x)lapply(1:length(o),function(y){if(o[y]==x)Landscapelist[y]})))
  new("LandscapeHistory",lo)
}

nbpar <- function(x) {unlist(lapply(x,function(x) switch(x[],
                            constant=1,
                            proportional=1,
                            enveloppe=2,
                            envelin=2,
                            quadratic=4,
                            fat_tail1=2,
                            gaussian=1,
                            exponential=1,
                            contiguous=1,
                            contiguous8=1,
                            island=1,
                            fat_tail2=2,
                            contiguous_long_dist_mixt=2,
                            gaussian_long_dist_mixt=2)
                            ))
                      }


validityNicheModel = function(object){
  if(class(object@variables)!="character")stop("error in NicheModel variables : variables just accept character!")
  if(!is.list(object@parameterList))stop("error in NicheModel parameterList : parameterList just accept list!")
  if(class(object@reactNorms)!="character")stop("error in NicheModel reactNorms : reactNorms just accept character!")
  if(FALSE%in%lapply(object@parameterList,is.numeric))stop("error in NicheModel parameter list : Parameter list just accept numeric!")
  if(length(object@variables)!=length(object@reactNorms))stop("error in NicheModel : number of variables and number of reaction norms do not correspond")
  notMatching <- (unlist(lapply(1:length(object@parameterList),function(x) nbpar(object@reactNorms[x]) != length(object@parameterList[[x]])))) 
  if (any(notMatching)) stop(paste("error in NicheModel : number of paremeters and reactionNorm do not match for variable ",which(notMatching),". ",sep=""))
  #              if grep("(",object@form)
  TRUE 
}

setClass("NicheModel",
         slots = c(variables="ANY",parameterList="ANY",reactNorms="ANY"),#,form="character"),
         validity=validityNicheModel
)

#setClass("form",
#		slots=c(numberVar="integer",additions="integer",multiplications="integer",openParenthesis="integer",closeParenthesis="integer"),
#      validity=function(object){
#        if (length(object@openPatrenthesis)!=length(object@closeParenthesis)) stop("error in 'form': different number of openParenthesis and closeparenthesis")
#        if (numberVar<=length(additions)+length(multiplications)) stop("error in 'form': too much operations")
#        if (order(append(addition,multiplications))) stop("error in 'form': too much operations")
#  }
#)


NicheModel<-function(variables=characterVector1,parameterList=listOfNumeric,reactNorms=characterVector2){#,form=formul){
  names(parameterList)=variables
  names(reactNorms)=variables
  new("NicheModel",variables=variables,parameterList=parameterList,reactNorms=reactNorms)#,form=form)
}

listOfMigrationShape<-c("fat_tail1","gaussian","exponential","contiguous","contiguous8","island","fat_tail2","contiguous_long_dist_mixt","gaussian_long_dist_mixt")

validityMigrationModel=function(object){
  if(!is.character(object@shapeDisp))stop("error in  MigrationModel shapeDisp : ShapeDisp just accept character!")
  if(length(object@shapeDisp)!=1)stop("error in  MigrationModel shapeDisp : ShapeDisp is  not valid because it contains more or less than one shape!")
  if(!object@shapeDisp%in%listOfMigrationShape)stop("error in  MigrationModel shapeDisp : the given shape not exist for MigrationModel !")
  if(FALSE%in%lapply(object@pDisp,is.numeric))stop("error in MigrationModel pDisp : pDisp just accept numeric!")
  if(nbpar(object@shapeDisp)!=length(object@pDisp))stop("error in MigrationModel : number of paremeters and shapeDisp do not match for variable")
  TRUE
}

setClass("MigrationModel",
         slots = c(shapeDisp="ANY",pDisp="ANY"),
         validity = validityMigrationModel
)



MigrationModel<-function(shape=character,param=p){
  new("MigrationModel",shapeDisp=shape,pDisp=param)
}

setClass("EnvDinModel",
         slots = c(K="ANY",R="ANY",migration="ANY"),
         validity=function(object){
           if(class(object@K)!="NicheModel")stop("Error in envDinModel K : K just accept NicheModel !")
           if(class(object@R)!="NicheModel")stop("Error in envDinModel R : R just accept NicheModel !")
           if(class(object@migration)!="MigrationModel")stop("Error in envDinModel migration : migration just accept MigrationModel !")
         }
         )

EnvDinModel<-function(K=nichemodelK,R=nichemodelR,migration=m){
  new("EnvDinModel",K=K,R=R,migration=migration)
}

setClass("TransitionBackward",
                             contains = "matrix",
                             validity = function(object){
                               if (all(nrow(object)==0))stop("The matrix is empty.")
                               if (nrow(object)!=ncol(object))stop("The matrix is not square")
                               if (!all(rowSums(object)>0.999999999 && rowSums(object)<1.111111111))stop("The sum of probabilities in each row is not 1.")
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

######### SET METHODS ##########################################################################################################

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
                   									       constant={setValues(object[[x]],rep(model@parameterList[[x]],ncell(object[[x]])))},
                   									       #proportional = {values(object[[x]])=object[[x]]*model@parameterList[[x]]},
                   									       enveloppe = {object[[x]]=enveloppe(object[[x]],model@parameterList[[x]])},
                   									       envelin={object[[x]]=envelinear(object[[x]],model@parameterList[[x]])},
                   									       conQuadratic={object[[x]]=conQuadratic(object[[x]],model@parameterList[[x]])}, 
                                                 #conquadraticskewed=conquadraticskewed(object[,,(model@variables==x)],p),
                                                 #conquadraticsq=conquadraticsq(object[,,(model@variables==x)],p),
                                                 #conquadraticskewedsq=conquadraticskewedsq(object[,,(model@variables==x)],p)
                   									       stop("This variable does not exist for NicheModel !")
                   									)
                         }
                )
            Y=prod(stack(Y))
          }
)



enveloppe <- function(X,p){
  if(length(p)!=2)stop("The parameter is  not valid because it contains more or less than two values")
  else X>=p[1]&X<=p[2]
}

envelinear <- function(X, p) {
  if(length(p)!=2)stop("The parameter is  not valid because it contains more or less than two values")
  else (X-p[1])/(p[2]-p[1])*enveloppe(X,p)
}
#plot(1:100,envelinear(1:100,c(20,60)))

constant <- function(X,p){X[]<-p}

conQuadratic <- function(X,p)
{
  if(length(p)!=2)stop("The parameter is  not valid because it contains more or less than two values")
  else -4*(X-p[2])*(X-p[1])/((p[2]-p[1])^2)*enveloppe(X,p)
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
######### CREER TRANSITION MATRIX ###############################################################################



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
    #if model["shapeDisp"]=="contiguous" matrix()
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
        return(migration)
  }
)

############### TEST DES FONCTION 
######### Landscape############
r1<- raster(ncol=2, nrow=1)
r1[] <- rep(1,2:2)
s <- stack(x=c(r1))
r0<- raster(ncol=0, nrow=0)
s <- stack(x=c(r0))
s<-1
r2<- raster(ncol=2, nrow=1)
r2[] <- rep(1,2:2)
s<- stack(x=c(r1,r2))

p<-c(as.Date("2000-01-10"),as.Date("2000-01-11"))
p<-"2000-01-11"
p<-as.Date("2000-01-11")
p<-c(as.Date("2000-01-11"),as.Date("2000-01-11"))
p<-c(as.Date("2000-01-11"),as.Date("2000-01-10"))
p<-c(as.Date("2000-01-11"),as.Date("2000-01-11"),as.Date("2000-01-11"))

v<-"a"
v<-c("a","b")
v<-c("layer.1","layer.2")
v<-1

lands<-Landscape(rasterstack = s,period = p,vars = v)

######### LandscapeHistory ################
r1<- raster(ncol=2, nrow=1)
r1[] <- rep(1,2:2)
r2<- raster(ncol=2, nrow=1)
r2[] <- rep(2,2:2)
s<- stack(x=c(r1,r2))
p1<-c(as.Date("2000-01-11"),as.Date("2000-01-12"))
p2<-c(as.Date("2000-01-13"),as.Date("2000-01-13"))
lands<-Landscape(rasterstack = s,period = p1,vars = v)
lands2<-Landscape(rasterstack = s,period = p2,vars = v)
lLands<-list(lands,lands2)
lLands<-c(lands2,lands)
lLands<-c(1,lands)
lLands<-c(lands,1)
lLands<-list(lands,lands)
landsH<-LandscapeHistory(lLands)
landsH[[1]]["period"]
#
######### NicheModel ###################
vari<-c("l")
vari<-"l"
vari<-c(1)
vari<-c("l",1)
vari<-c(1,"l")
vari<-c("l","t")
vari<-c(NA)
vari<-c(NULL)

para<-list(1)
para<-list(c(1,5))
para<-list(1,c(1,5))
para<-list(1,2)
para<-list(c(2,3),c(1,5))
para<-list(1,c(1,5))

rea<-c(l="constant")
rea<-c(l="constant",t="envelin")
rea<-c(l="constant",t="enveline")
rea<-c(l="constant",t="envelin",p="enveloppe")

model<-NicheModel(variables=vari,parameterList=para,reactNorms=rea)
model



######### MigrationModel ###################
migraShape<-"fat_tail1"
migraShape<-"gaussian"
migraShape<-"exponential"
migraShape<-"contiguous"
migraShape<-"contiguous8"
migraShape<-"island"
migraShape<-"fat_tail2"
migraShape<-"contiguous_long_dist_mixt"
migraShape<-"gaussian_long_dist_mixt"
migraShape<-"constant"
migraShape<-"yo"
migraShape<-c("fat_tail1","gaussian")
migraShape<-1

migrapDisp<-1
migrapDisp<-c(1)
migrapDisp<-c(1,2)
migrapDisp<-"1"
migrapDisp<-1.2
migrapDisp<-c(1.2)
migrapDisp<-c(1,"1")


migramodel<-MigrationModel(shape = migraShape,param = migrapDisp)


######### EnvDinModel ###################
vari<-c("l","t")

para<-list(1,c(1,5))
para<-list(1,2)
para<-list(c(2,3),c(1,5))

rea<-c(l="constant",t="envelin")

migrapDisp<-c(1,2)
migraShape<-"fat_tail1"

nicheK<-NicheModel(variables = vari, parameterList = para,reactNorms = rea)
nicheK<-1

nicheR<-NicheModel(variables = vari, parameterList = para,reactNorms = rea)
nicheR<-"2"
  
migramodel<-MigrationModel(shape = migraShape,param = migrapDisp)
migramodel<-NicheModel(variables = vari, parameterList = para,reactNorms = rea)
env1<-EnvDinModel(K=nicheK,R=nicheR,migration = migramodel)



######### test fonction ###################
r1<- raster(ncol=2, nrow=1)
r1[] <- rep(1,2:2)
r2<- raster(ncol=2, nrow=1)
r2[] <- rep(1,2:2)
s<- stack(x=c(r1,r2))

p1<-as.Date("2000-01-11")
v<-c("l","p")
lands<-Landscape(rasterstack = s,period = p1,vars = v)
lands<-Landscape(rasterstack = s2,period = p1,vars = v)


######### CONSTRUCTION D'UN TRANSITION MATRIX ################
r1<- raster(ncol=2, nrow=2)
r1[] <- rep(2:5,1)
r2<- raster(ncol=2, nrow=2)
r2[] <- rep(2,2:2)
s<- stack(x=c(r1,r2))
p1<-as.Date("2000-01-11")
variK<-c("l","t")
variR<-c("c","t")
paraK<-list(c(0,5),2)
paraR<-list(2,2)
reaK<-c(l="envelin",t="constant")
reaR<-c(l="constant",t="constant")
extent(s)<-c(0,2,0,2)
lscp1<-Landscape(rasterstack = s,period=p1,vars=vari)
modelK<-NicheModel(variables=vari,parameterList=para,reactNorms=rea)
modelR<-NicheModel(variables=vari,parameterList=para,reactNorms=rea)
m<-MigrationModel(shape="gaussian",param = (1/1.96))
a<-runNicheModel(lscp1,modelK)
migrationMatrix(lscp1,m)

edm1<-EnvDinModel(K=modelK,R=modelR,migration = m)
runEnvDinModel(lscp1,edm1)
transi1<-createTransitionMatrix(lscp1,edm1)
