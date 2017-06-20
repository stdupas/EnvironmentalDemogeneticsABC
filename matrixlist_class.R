########MANIPULATION CLASS####################################################################
mat3 <- rbind(c(0.2,0.8),c(0.6,0.4))
mat5 <- rbind(c(0.2,0.8),c(0.1,0.9))
d <- as.Date( c("2017-02-01","2017-02-01"))
d1<- as.Date( c("2017-02-01","2017-09-01"))
v <- "poid"
v1<-"taille"
wa<-LandPeriod(matrix=mat3,date=d,var=v)
wo<-LandPeriod(matrix=mat5,date=d1,var=v)
l=list(wa,wo)
l
ya=new("MatrixList",Liste=l)
yo=MatrixList(wa)
lengthM(ya)

dimnames(mat3)
dim(mat3)
names(mat3)
mLength(ya)
dimnamesM(ya)
dimM(ya)
namesM(ya)

############class#############################################################################
LandPeriod<-setClass("LandPeriod",
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

MatrixList<-setClass("MatrixList",
                     contains = "list",
                     slots = c(Liste="list"),
                     validity = function(object){
                       if(class(object@Liste[[1]])[1]!="LandPeriod")stop("The element is not a LandPeriod")
                     }
                     )


MatrixList<-function(mat=mat3){
  Lis=list(mat)
  new("MatrixList",Liste=Lis)
}

############METHOS ~ MATRIX###########################################################################
setGeneric(
  name = "getList",
  def=function(object){return(standardGeneric("getList"))}
)

setMethod(
  f="getList",
  signature = "MatrixList",
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
  signature = "MatrixList",
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
  signature = "MatrixList",
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
  signature = "MatrixList",
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
  signature = "MatrixList",
  definition = function(object){
    lapply(lapply(getList(object),getMatrix),names)
  }
)

############METHODS###################################################################################

setGeneric(
  name = "getMatrix",
  def=function(object){return(standardGeneric("getMatrix"))}
)

setMethod(
  f="getMatrix",
  signature = "LandPeriod",
  definition = function(object){
    return(object@matrix)
  }
)

setGeneric(
  name = "getDate",
  def=function(object){return(standardGeneric("getDate"))}
)

setMethod(
  f="getDate",
  signature = "LandPeriod",
  definition = function(object){
    return(object@date)
  }
)


setGeneric(
  name = "getVar",
  def=function(object){return(standardGeneric("getVar"))}
)

setMethod(
  f="getVar",
  signature = "LandPeriod",
  definition = function(object){
    return(object@var)
  }
)


setGeneric(
  name = "setMatrix",
  def=function(object,mat){return(standardGeneric("setMatrix"))}
)

setGeneric(
  name = "setMatrixb<-",
  def=function(object,value){return(standardGeneric("setMatrixb<-"))}
)

setReplaceMethod(
  f = "setMatrixb" ,
  signature = "LandPeriod" ,
  definition = function ( object , value ){
      object@matrix <- value
    }
)

setMethod(
  f="setMatrix",
  signature = "LandPeriod",
  definition = function(object,mat){
    (`setMatrixb<-`(object,mat))
  }
)

setGeneric(
  name = "setDate",
  def=function(object,date){return(standardGeneric("setDate"))}
)

setGeneric(
  name = "setDateb<-",
  def=function(object,value){return(standardGeneric("setDateb<-"))}
)

setReplaceMethod(
  f = "setDateb" ,
  signature = "LandPeriod" ,
  definition = function ( object , value ){
    object@date <- value
  }
)

setMethod(
  f="setDate",
  signature = "LandPeriod",
  definition = function(object,date){
    `setDateb<-`(object,date)
  }
)


setGeneric(
  name = "setVar<-",
  def=function(object,value){return(standardGeneric("setVar<-"))}
)

setReplaceMethod(
  f = "setVar" ,
  signature = "LandPeriod" ,
  definition = function (object, value ){
    assign(object@var,value,envir = e1)
    return(invisible())
  }
)

