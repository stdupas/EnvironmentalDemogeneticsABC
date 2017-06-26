########MANIPULATION CLASS####################################################################
mat3 <- rbind(c(0.2,0.8),c(0.6,0.4))
mat5 <- rbind(c(0.2,0.8),c(0.1,0.9))
d <- as.Date( c("2017-02-01","2017-02-01"))
d1<- as.Date( c("2017-02-01","2017-09-01"))
v <- "poid"
v1<-"taille"
wa<-LandPeriod(matrix=mat3,date=d,var=v)
wo<-LandPeriod(matrix=mat5,date=d1,var=v)
wa["var"]<-v1
wa["var"]
wa
l=list(wa,wo)
l
ya<-LandPeriodList(wa)
l1=list(wo)
ya["Liste"]<-l1
ya

l1yo=LandPeriodList(wa)
lengthM(ya)
addMatrix(ya,wo)
dimnames(mat3)
dim(mat3)
names(mat3)
mLength(ya)
dimnamesM(ya)
dimM(ya)
namesM(ya)
ya
setMatrixb(wa,mat5)
`setMatrix<-`(wa,mat5)
class(wa)
setVar(wa)<-v
wa
`setDateb<-`()

############class#############################################################################
LandPeriod<-setClass("LandPeriod",
                     #contains = "matrix",
                     slots = c(matrix= "matrix", date="Date",var="character"),
                     validity = function(object){
                       if(is.matrix(object@matrix)){}else{stop("The matrix is not a matrix")}
                       if(length(object@date)>2)stop("To many date ")
                       if(object@date[[1]]>object@date[[length(object@date)]])stop("The begin is after the end")
                       if(is.character(object@var)){}else{stop("var is not character")}
                     }
)

Landperiod<-function(mat=mat3,da=d,va=v){
  new("LandPeriod",matrix=ma,date=da,var=va)
}

LandPeriodList<-setClass("LandPeriodList",
                     contains = "list",
                     #slots = c(Liste="list"),
                     validity = function(object){
                       if(class(object@Liste[[1]])[1]!="LandPeriod")stop("The element is not a LandPeriod")
                     }
                     )


LandPeriodList<-function(mat=mat3){
  Lis=list(mat)
  new("LandPeriodList",Liste=Lis)
}

############METHOS ~ MATRIX###########################################################################

setGeneric(
  name = "addLandPeriod",
  def=function(object,Lp){return(standardGeneric("addLandPeriod"))}
)

setMethod(
  f="addLandPeriod",
  signature = "LandPeriodList",
  definition = function(object,Lp){
    object[]<-Lp
    return(object[])
  }
)


setGeneric(
  name = "getList",
  def=function(object){return(standardGeneric("getList"))}
)

setMethod(
  f="getList",
  signature = "LandPeriodList",
  definition = function(object){
    return(object@Liste)
  }
)

setGeneric(
  name = "lengthM",
  def=function(object){return(standardGeneric("lengthM"))}
)

setMethod(
  f="lengthM",
  signature = "LandPeriodList",
  definition = function(object){
    lapply(lapply(getList(object),getMatrix),length)
  }
)

setGeneric(
  name = "dimnamesM",
  def=function(object){return(standardGeneric("dimnamesM"))}
)

setMethod(
  f="dimnamesM",
  signature = "LandPeriodList",
  definition = function(object){
    lapply(lapply(getList(object),getMatrix),dimnames)
  }
)
        
setGeneric(
  name = "dimM",
  def=function(object){return(standardGeneric("dimM"))}
)

setMethod(
  f="dimM",
  signature = "LandPeriodList",
  definition = function(object){
    lapply(lapply(getList(object),getMatrix),dim)
  }
)


setGeneric(
  name = "namesM",
  def=function(object){return(standardGeneric("namesM"))}
)

setMethod(
  f="namesM",
  signature = "LandPeriodList",
  definition = function(object){
    lapply(lapply(getList(object),getMatrix),names)
  }
)

############METHODS###################################################################################
setMethod(
  f ="[",
  signature = c(x="LandPeriod" ,i="character",j="missing"),
  definition = function (x ,i ,j , drop ){
    switch ( EXPR =i,
             "matrix" ={return(x@matrix)} ,
             "date" ={return(x@date)} ,
             "var" ={return(x@var)} ,
             stop("This slots doesn't exist!")
    )
  }
)


setReplaceMethod(
  f ='[',
  signature = c(x="LandPeriod" ,i="character",j="missing",value="ANY") ,
  definition = function(x,i,j,value){
    switch(EXPR=i,
           "matrix" ={x@matrix<-value} ,
           "date" ={x@date<-value} ,
           "var" ={x@var<-value} ,
           stop("This slots doesn't exist!")
    )
    validObject(x)
    return(x)
  }
)

setMethod(
  f ="[",
  signature = c(x="LandPeriodList" ,i="missing",j="ANY"),
  definition = function (x ,i ,j , drop )return(x@Liste)
)


setReplaceMethod(
  f ='[',
  signature = c(x="LandPeriod" ,i="character",j="missing",value="list") ,
  definition = function(x,i,j,value){
    switch(EXPR=i,
           "Liste" ={x@Liste<-value} ,
           stop("This slots doesn't exist!")
    )
    validObject(x)
    return(x)
  }
)



setReplaceMethod(
  f ='[',
  signature = c(x="LandPeriodList" ,i="missing",j="missing",value="list") ,
  definition = function(x,i,j,value){ 
    x@Liste<-value
    validObject(x)
    return(x)
  }
)






