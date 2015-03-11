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
  nCell <- length(latitude)*length(longitude)/(aggregationParam^2)
  EnvData <- array(NA,
                   dim=c(nCell,length(Dates),length(ncFiles)),
                   dimnames=list(1:nCell,as.character(Dates),1:length(ncFiles)))
  pixels <- SpatialPixels(SpatialPoints(cbind(longitude,latitude)))
  rasterLayer <- raster(pixels)
  rasterAgg <-  aggregate(rasterLayer,aggregationParam)
  rasterStackAgg <- stack(rasterAgg)
  layers <- NA
  for (ncFile in ncFiles)
  {
    variable = dimnames(EnvData)[[3]][which(ncFiles==ncFile)] = layers[which(ncFiles==ncFile)] = strsplit(ncFile,"_")[[1]][1]
    nc=open.nc(paste(ncDirectory,ncFile,sep=""))
    enVarNC <- var.get.nc(nc,variable)
    for (i in 1:length(Dates))
    {
      values(rasterLayer) <- matrix(enVarNC[,,i],nrow=length(latitude),ncol=length(longitude),byrow=TRUE)
      rasterAgg <- aggregate(rasterLayer,aggregationParam)
      EnvData[,as.character(Dates[i]),variable] <- values(rasterAgg)
    }
    rasterStackAgg <- stack(rasterStackAgg,rasterAgg)
  }
  names(rasterStackAgg) <- layers
  list(EnvData,rasterStackAgg)
}
