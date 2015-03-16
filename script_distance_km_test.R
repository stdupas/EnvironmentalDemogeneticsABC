library("argosfilter")

distanceMatrixFromRaster =
function(object){
        # Computes a pairwise distance matrix from a raster object
        #
        # Args:
        #   object: a raster object from which computes distances
        #
        # Returns:
        #   A matrix of distances in meters if a coordinate system is precised
        
        # Extract coordinates from raster object
        coords = xyFromCell(object = object, cell = 1:ncell(object), spatial=FALSE)
        lat = coords[,1]
        lon = coords[,2]
        # Compute distance matrix
        dist = NULL
        for(i in 1:length(lat)) {
            res = NULL
            for(j in 1:length(lon)) {
                res = c(res,distance(lat[i],lat[j],lon[i],lon[j]))
            }
            dist = rbind(dist,res)
        }
        return(dist)
}