setMethod(
  f="getMatrix",
  signature = "TransitionBackward",
  definition = function(object){return(object@matrix)}
)
existsMethod(f="getMatrix",signature = "TransitionBackward")

setMethod(
  f="setMatrix",
  signature = "TransitionBackward",
  definition = function(object,matrix){
    object@matrix<-matrix
    return(object@matrix)}
)

setMethod(
  f = "laplaceMatrix",
  signature = "TransitionBackward",
  definition = function(object){
    matrixD = diag(rep(1,dim(object@matrix)[1])) # diagonal equals to 1
    laplacianMatrix = matrixD - object@matrix
    laplacianMatrix[is.na(laplacianMatrix)]<-0 # replace NA by 0
    cat("laplacian",laplacianMatrix)
  }
)

existsMethod(f="laplaceMatrix",signature = "TransitionBackward")



setMethod(
  f="ordinary_laplacian",
  signature = "TransitionBackward",
  definition = function(object){
    markovB<-new("markovchain", states=dimnames(transition)[[1]], transitionMatrix=transition)
    PI<-diag(steadyStates(markovB)[1,])
    PI - PI%*%transition
  }
)



setMethod(
  f="commute_time_undigraph ",
  signature = "TransitionBackward",
  definition = function(object){
    laplacian = laplaceMatrix(object)
    inverseMP = ginv(laplacian) # generalized inverse matrix  (Moore Penrose)
    diag = diag(inverseMP) # get diagonal of the inverse matrix
    mii = matrix(diag, nrow =dim(inverseMP), ncol = dim(inverseMP))
    mjj = t(mii)
    mij = inverseMP
    mji = t(mij)
    commute_time = mii + mjj - mij - mji
    commute_time
    }
)



setMethod(
  f="hitting_time_digraph",
  signature = "TransitionBackward",
  definition = function(object){
    Ones <- rep(1,dim(getMatrix(object))[1])
    markovB<-new("markovchain", states=dimnames(getMatrix(object))[[1]], transitionMatrix=getMatrix(object))
    pi_<-steadyStates(markovB)[1,]
    PI <- diag(pi_)
    L <- PI - PI%*%getMatrix(object)
    Z <- ginv(L + pi_%*%t(pi_))
    H <- Ones%*%t(diag(Z))-Z
    H
  }
)

setMethod(
  f="genetDistUndigraph",
  signature = "TransitionBackward",
  definition = function(object,popSize,mutation_rate){
    commute_time <- commute_time_undigraph(object)
    #genetic_dist = commute_time / (8* popSize)
    #genetic_dist = commute_time / (8* (sum(popSize)/(dim(popSize)[1]*dim(popSize)[2])))
    genetDist =  (commute_time/4+2*sum(valuesA(popSize))) * mutation_rate
    genetDist
  }
)


setMethod(
  f="genetDistDigraph",
  signature = "TransitionBackward",
  definition = function(object,popSize,mutation_rate,method="Goldstein95"){
    H <- hitting_time_digraph(object)
    dim2 <- dim(H);dim2[[3]]=2
    H2 <- array(c(H,t(H)),dim=dim2)
    MinH <- apply(H2,c(1,2),min)
    #H2[,,1] <- MinH; H2[,,2] = (H+t(H))/2
    #MinH2 <- apply(H2,c(1,2),min)
    genetic_dist = (MinH/2 + 2*sum(valuesA(popSize)))* mutation_rate
    genetic_dist
  }
)


############################ METHODS NON ADAPTÉE #############################################

setMethod(
  f="degree2km",
  signature = "",
  definition = function(rasterStack){
    x_origin = ((xmin(rasterStack)+xmax(rasterStack))/2) #longitude origin
    y_origin = ((ymin(rasterStack)+ymax(rasterStack))/2) #latitude origin
    x_destination = (x_origin + xres(rasterStack)) #longitude of destination point
    y_destination = (y_origin + yres(rasterStack)) #latitude of destination point
    
    dist_degree <- acos(sin(x_origin)*sin(x_destination)+cos(x_origin)*cos(x_destination)*cos(y_origin-y_destination))
    dist_km = dist_degree * 111.32
    dist_km
  }
)


setMethod(
  f="Aggregate_and_adjust_raster_to_data ",
  signature = "",
  definition = function(Envir_raster_stack,release,recovery,extend_band_size,aggregate_index){
    samples <- SpatialPoints(rbind(na.omit(release[,c("X","Y")]),na.omit(recovery[,c("X","Y")])))
    if (aggregate_index > 1) {
      Envir_raster_stack <- aggregate(crop(Envir_raster_stack,extent(samples)+extend_band_size), fact=aggregate_index, fun=mean, expand=TRUE, na.rm=TRUE)
      } 
    else {
      Envir_raster_stack <- crop(Envir_raster_stack,extent(samples)+extend_band_size)
      }
    Envir_raster_stack
  }
)

setMethod(
  f="conquadraticskewed",
  signature = "",
  definition = function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12")))){
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


setMethod(
  f="conquadraticskewedsq",
  signature = "",
  definition = function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12")))){
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
)

setMethod(
  f="conquadratic",
  signature = "",
  definition = function(X,p){
    xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
    xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]  
    yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]  
    (yopt-(4*yopt/(xmax-xmin)^2)*(X-(xmin+xmax)/2)^2)*((X>xmin)&(X<xmax))
  }
)


setMethod(
  f="conquadraticsq",
  signature = "",
  definition = function(X,p){
    xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
    xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]  
    yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]  
    res = (yopt-(4*yopt/(xmax-xmin)^2)*(X-(xmin+xmax)/2)^2)*((X>xmin)&(X<xmax))
    res + res * (1-res)
  }
)


setMethod(
  f="enveloppe",
  signature = "",
  definition = function(X,p){
    p[rep("Yopt",dim(X)[1]),colnames(X)]*((X>p[rep("Xmin",dim(X)[1]),])&(X<p[rep("Xmax",dim(X)[1]),colnames(X)]))
  }
)

setMethod(
  f="envelinear",
  signature = "",
  definition = function(X,p,log=FALSE){
    Yxmin = p[rep("Yxmax",dim(X)[1]),colnames(X)]
    Yxmax = p[rep("Yxmin",dim(X)[1]),colnames(X)]
    Xmin=  p[rep("Xmin",dim(X)[1]),colnames(X)]
    Xmax=  p[rep("Xmax",dim(X)[1]),colnames(X)]
    a = (Yxmin - Yxmax) / (Xmin - Xmax)
    b = Yxmin - Xmin * a
    if (log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
  }
)



setMethod(
  f="envelin0",
  signature = "",
  definition = function(X,p,log=FALSE){
    Yxmin = p[rep("Yxmax",dim(X)[1]),colnames(X)]
    Yxmax = p[rep("Yxmin",dim(X)[1]),colnames(X)]
    Xmin=  p[rep("Xmin",dim(X)[1]),colnames(X)]
    Xmax=  p[rep("Xmax",dim(X)[1]),colnames(X)]
    a = (Yxmin - Yxmax) / (Xmin - Xmax)
    b = Yxmin - Xmin * a
    if (log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
  }
)


setMethod(
  f="proportional",
  signature = "",
  definition = function(X,p,Log=FALSE){
    if (Log) {log(p[rep("a",dim(X)[1]),colnames(X)]*X)} else {
      p[rep("a",dim(X)[1]),colnames(X)]*X
    }
  }
)


setMethod(
  f="",
  signature = "",
  definition = function(){
    
  }
)