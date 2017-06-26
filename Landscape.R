##########MANIPULATION CLASS################################


ar<-array(1:12,c(2,2,3))
per<-as.Date( c("2017-02-01","2017-02-01"))
vari<-c("l","t","p")
para<-list(c(1,3),2,c(4,5.2))
rea<-c(l="constant",t="constant",p="constant")
dim(lscp)[3]
lscp<-Landscape(array=ar,period=per,listvar=vari)
lscp
length(lscp@listVar)

is.Landscape(lscp)
lista<-list(lscp)
lh1<-LandscapeHistory(lista)
a=list(as.numeric(1:2))
b=c(a,a)
class(a[1])
lapply(b,is.numeric)


p=c("Xmin"=10,"Xmax"=20,"Xopt"=18,"Ymax"=0.1)
p[rep("a",dim(ar[,,3])[1]),colnames(ar[,,3])]*ar[,,3]

dim(ar[,,3])


model<-NicheModelCombinedPerPeriod(vari,para,rea,per)
model

nicheModelVar(lscp,model)
lscp

##########SET CLASS#########################################
Landscape<-setClass("Landscape",
         contains = "array",
         slots = c(period="Date",listVar="character"),
         validity = function(object){
           if(object@period[1]>object@period[2])stop("A period begin after it end")
           if(length(object@listVar)!=dim(object)[3])stop("Not have 1 variable for each layer")
         }
         )


Landscape<-function(array=a,period=p,listvar=l){
  if(length(period)>1){
    for(i in 1:length(period)/2){
      if(period[1]!=period[2*i-1]&&period[2]!=period[2*i])stop("The period is not the same for each layer")
    }
    peri<-c(period[1:2])
  }else{peri<-c(period[1],period[1])}
  new("Landscape",array,period=peri,listVar=listvar)
}


LandscapeHistory<-setClass("LandscapeHistory",
                           contains = "list",
                           validity = function(object){
                             if(is.Landscape(object)){}
                             else stop("An element of the list is not a Landscape")
                           }
                           )
LandscapeHistory<-function(a){
  new("LandscapeHistory",a)
}


NicheModelCombinedPerPeriod<-setClass("NicheModelCombinedPerPeriod",
                                      slots = c(variable="character",parameter="list",reactNormList="character",period="Date"),
                                      validity = function(object){
                                        if(FALSE%in%lapply(object@parameter,is.numeric))stop("Parameter just accept numeric!")
                                        if(object@period[1]>object@period[length(object@period)])stop("A period begin after it end")
                                        if(length(object@variable)!=length(object@reactNormList))stop("Each variable don t have 1 reactNorm")
                                      }
)

NicheModelCombinedPerPeriod<-function(variable=var,parameter=par,reactNormList=rea,period=pe){
  new("NicheModelCombinedPerPeriod",variable=variable,parameter=parameter,reactNormList=reactNormList,period=period)
}

#######SET METHODS##########################################

setGeneric(
  name = "validDate",
  def=function(object){return(standardGeneric("validDate"))}
)

setMethod(
  f="validDate",
  signature = "Landscape",
  definition = function(object){
    dime<-c(1:(length(object@period)/2))
    j=1
    for (i in dime){
      if(object@period[j+1]<object@period[j]){return(FALSE)}
      j=j+2
    }
    return(TRUE)
  }
)


setGeneric(
  name = "is.Landscape",
  def=function(object){return(standardGeneric("is.Landscape"))}
)

setMethod(
  f="is.Landscape",
  signature = "Landscape",
  definition = function(object){
    if(class(object)[1]=="Landscape")return(TRUE)
    else return(FALSE)
  }
)

setMethod(
  f="is.Landscape",
  signature = "LandscapeHistory",
  definition = function(object){
    if(FALSE%in%lapply(object,is.Landscape)){return(FALSE)}
    else{return(TRUE)}
  }
)



setGeneric(
  name = "nicheModelVar",
  def=function(object,model){return(standardGeneric("nicheModelVar"))}
)



setMethod("nicheModelVar",
          signature=c("Landscape","NicheModelCombinedPerPeriod"),
          definition = function(object,model){                  #X=object, p=,shape=
            lapply(model@variable,function(x) switch(model@reactNormList[[x]],
                                                 constant=model@parameter[model@variable==x],
                                                 proportional = proportional(object[,,(model@variable==x)],p), 
                                                 linear = linear(object[,,(model@variable==x)],p),
                                                 enveloppe=enveloppe(object[,,(model@variable==x)],p),
                                                 envelin=envelinear(object[,,(model@variable==x)],p),
                                                 envloglin=envelinear(object[,,(model@variable==x)],p,log=TRUE),
                                                 loG = log(object[,,(model@variable==x)]),
                                                 conquadratic=conquadratic(object[,,(model@variable==x)],p),
                                                 conquadraticskewed=conquadraticskewed(object[,,(model@variable==x)],p),
                                                 conquadraticsq=conquadraticsq(object[,,(model@variable==x)],p),
                                                 conquadraticskewedsq=conquadraticskewedsq(object[,,(model@variable==x)],p)
                           )
                )
          }
)


setGeneric(
  name = "proportional",
  def=function(X,p,Log=FALSE){return(standardGeneric("proportional"))}
)


setMethod(
  f="proportional",
  signature=c("matrix"),
  definition=function(X,p){
    if(Log){log(p[ rep ("a",dim(X)[1]),colnames(X)]*X)}else {
      p[rep("a",dim(X)[1]),colnames(X)]*X
    }
  }
)


setGeneric(
  name = "enveloppe",
  def=function(X,p){return(standardGeneric("enveloppe"))}
)


setMethod(
  f="enveloppe",
  signature="matrix",
  definition=function(X,p){
    p[rep("Yopt",dim(X)[1]),colnames(X)]*((X>p[rep("Xmin",dim(X)[1]),])&(X<p[rep("Xmax",dim(X)[1]),colnames(X)]))
  }
)  



setGeneric(
  name = "conquadraticsq",
  def=function(X,p){return(standardGeneric("conquadraticsq"))}
)



setMethod(
  f="conquadraticsq",
  signature="",
  definition=function(X,p)
  {
    xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
    xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]  
    yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]  
    res = (yopt-(4*yopt/(xmax-xmin)^2)*(X-(xmin+xmax)/2)^2)*((X>xmin)&(X<xmax))
    res + res * (1-res)
  }
)

setGeneric(
  name = "conquadraticskewedsq",
  def=function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12")))){return(standardGeneric("conquadraticskewedsq"))}
)


setMethod(
  f="conquadraticskewedsq",
  signature="",
  definition=function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))))
  {
    Yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]
    Xopt = p[rep("Xopt",dim(X)[1]),colnames(X)]
    Xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]
    Xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
    alpha <- -log(2)/log((Xopt-Xmin)/(Xmax-Xmin))
    Xprime<- ((X-Xmin)/(Xmax-Xmin))^alpha*(Xmax-Xmin)+Xmin
    y <- (Yopt-4*Yopt/((Xmax-Xmin)^2)*(Xprime-(Xmin+Xmax)/2)^2)*((X>=Xmin)&(X<=Xmax))
    y+y*(Yopt-y)
  }
)


setGeneric(
  name = "conquadraticskewed",
  def=function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12")))){return(standardGeneric("conquadraticskewed"))}
)


setMethod(
  f="conquadraticskewed",
  signature="",
  definition=function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))))
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
  
)


setGeneric(
  name = "conquadratic",
  def=function(X,p){return(standardGeneric("conquadratic"))}
)

setMethod(
  f="conquadratic",
  signature="",
  definition=function(X,p)
  {
    xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
    xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]  
    yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]  
    (yopt-(4*yopt/(xmax-xmin)^2)*(X-(xmin+xmax)/2)^2)*((X>xmin)&(X<xmax))
  }
)


setGeneric(
  name = "envelinear",
  def=function(X,p,log=FALSE){return(standardGeneric("envelinear"))}
)

setMethod(
  f="envelinear",
  signature="",
  definition=function(X,p,log=FALSE)
  {
    Yxmin = p[rep("Yxmax",dim(X)[1]),colnames(X)]
    Yxmax = p[rep("Yxmin",dim(X)[1]),colnames(X)]
    Xmin=  p[rep("Xmin",dim(X)[1]),colnames(X)]
    Xmax=  p[rep("Xmax",dim(X)[1]),colnames(X)]
    a = (Yxmin - Yxmax) / (Xmin - Xmax)
    b = Yxmin - Xmin * a
    if (log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
  }
  
)



setGeneric(
  name = "linear",
  def=function(X,p,Log=FALSE){return(standardGeneric("linear"))}
)

setMethod(
  f="linear",
  signature="",
  definition=function(X,p,Log=FALSE)
  {
    Yx1 = p[rep("Yx1",dim(X)[1]),colnames(X)]
    Yx0 = p[rep("Yx0",dim(X)[1]),colnames(X)]
    a = (Yx1 - Yx0)
    b = Yx0
    if (Log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
  }
)  




















