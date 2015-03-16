nc2EnvDataAndRasterStack <- function(ncDirectory=paste(wd,"ForwardSimulData/",sep=""),aggregationParam=4)
{
  ncFiles=grep(".nc",list.files(ncDirectory),value=TRUE)
  nc=open.nc(paste(ncDirectory,ncFiles[1],sep=""))
  
  Dates <-   strsplit(ncFiles[1],"_")[[1]][length(strsplit(ncFiles[1],"_")[[1]])]
  Dates <- strsplit(strsplit(Dates,"[.]")[[1]][1],"-")[[1]]
  Dates <- data.frame(Annees = substr(Dates,1,4),mois = substr(Dates,5,6),jours = substr(Dates,7,8))
  Dates <- transform(Dates, Date = paste(Annees, mois, jours, sep='/'))[-c(1:3)]
  Dates <- as.Date(as.Date(as.character(Dates[1,1])):as.Date(as.character(Dates[2,1])),origin="1970-01-01")
  longitude <- var.get.nc(nc,'longitude')
  latitude <- var.get.nc(nc,'latitude')
  pixels <- SpatialPixels(SpatialPoints(cbind(longitude,latitude)))
  rasterLayer <- raster(pixels)
  if(aggregationParam>1) {
    if (any((dim(rasterLayer)%%aggregationParam)[1:2]!=c(0,0))) {
      cellsToRemove <- dim(rasterLayer)[1:2]%%aggregationParam 
      rasterLayer2 <- crop(rasterLayer,extent(as.vector(extent(rasterLayer))-c(0,res(rasterLayer)[1]*(cellsToRemove[1]),0,res(rasterLayer)[2]*(cellsToRemove[2]))))
    } else {rasterLayer2 <- rasterLayer}
    rasterAgg <- aggregate(rasterLayer2,aggregationParam)} else {rasterAgg <- rasterLayer}

  nCell <- ncell(rasterAgg) #length(latitude)*length(longitude)/(aggregationParam^2)
  EnvData <- array(NA,
                   dim=c(nCell,length(Dates),length(ncFiles)),
                   dimnames=list(1:nCell,as.character(Dates),1:length(ncFiles)))
  rasterStackAgg <- stack(rasterAgg)
  layers <- NA
  for (ncFile in ncFiles)#ncFile=ncFiles[1]
  {
    variable = dimnames(EnvData)[[3]][which(ncFiles==ncFile)] = layers[which(ncFiles==ncFile)] = strsplit(ncFile,"_")[[1]][1]
    nc=open.nc(paste(ncDirectory,ncFile,sep=""))
    enVarNC <- var.get.nc(nc,variable)
    cropBeforeAggregate <- any((dim(rasterLayer)%%aggregationParam)[1:2]!=c(0,0))
    for (i in 1:length(Dates))
    {
      values(rasterLayer) <- matrix(enVarNC[,,i],nrow=length(latitude),ncol=length(longitude),byrow=TRUE)
      if(aggregationParam>1) {
        if (cropBeforeAggregate) {
          rasterLayer2 <- crop(rasterLayer,extent(as.vector(extent(rasterLayer))-c(0,res(rasterLayer)[1]*(dim(rasterLayer)[1]%%aggregationParam),0,res(rasterLayer)[2]*(dim(rasterLayer)[2]%%aggregationParam))))
        }
        rasterAgg <- aggregate(rasterLayer2,aggregationParam)} else {rasterAgg <- rasterLayer}
      EnvData[,as.character(Dates[i]),variable] <- values(rasterAgg)
    }
    rasterStackAgg <- stack(rasterStackAgg,rasterAgg)
  }
  names(rasterStackAgg) <- layers
  list(EnvData,rasterStackAgg)
}
