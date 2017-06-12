setMethod(
  f = "ncellA",
  signature = "RasterLayer",
  definition = function(object){
    length(na.omit(values(object)))
  }
)
setMethod(
  f = "ncellA",
  signature = "RasterStack",
  definition = function(object){
    ncellA(object[[1]])
  }
)
setMethod(
  f = "ncellA",
  signature = "RasterBrick",
  definition = function(object){
    ncellA(object[[1]])
  }
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

setMethod(
  f = "cellNumA",
  signature = "RasterBrick",
  definition = function(object){
    cellNumA(object[[1]])
  }
)

setMethod(
  f = "xyFromCellA",
  signature = "RasterLayer",
  definition = function(object){
    xyFromCell(object,cellNumA(object))
  }
)
setMethod(
  f = "xyFromCellA",
  signature = "RasterStack",
  definition = function(object){
    xyFromCell(object,cellNumA(object))
  }
)
setMethod(
  f = "xyFromCellA",
  signature = "RasterBrick",
  definition = function(object){
    xyFromCell(object,cellNumA(object))
  }
)

setMethod(
  f = "valuesA",
  signature = "RasterLayer",
  definition = function(object){
    object[cellNumA(object)]
  }
)

setMethod(
  f = "valuesA",
  signature = "RasterStack",
  definition = function(object){
    values(object)[cellNumA(object),]
  }
)

setMethod(
  f = "valuesA",
  signature = "RasterBrick",
  definition = function(object){
    values(object)[cellNumA(object),]
  }
)

setMethod(
  f="get",
  signature = "spatialGenetic",
  definition = function(x) {
    cbind(x=x@x,y=x@y,Cell_numbers=x@Cell_numbers,data.frame(x))
  }
)


setMethod(
  f="coordinates",
  signature = "spatialGenetic",
  definition = function(obj) {
    cbind(x=obj@x,y=obj@y)
  }
)

setMethod(
  f="colnames",
  signature = "spatialGenetic",
  definition = function(x){
    colnames(get(x))
  })

setMethod(
  f="cellNumA",
  signature = "spatialGenetic",
  definition = function(object) {
    cbind(Cell_numbers=object@Cell_numbers)
  }
)



setMethod(
  f="aggregate",
  signature = "spatialGenetic",
  definition = function(x, by="Cell_numbers", FUN="mean") {
    cols <- colnames(x)
    X <- get(x)
    cells = as.integer(levels(as.factor(cellNumA(x))))
    Y <- X[1,];Y<-Y[-1,]
    for (cell in cells){
      for (j in cols){
        Y[cell,j] <- switch(FUN,
                            "mean"=mean(na.omit(X[X[by]==cell,j])))
      }      
    }
    spatialGenetic(Y[cells,])
  }
)
?dist
