CreateGenetArray <- function(rasK, nb_locus, initial_locus_value, Option, nind)
{
  # This function constructs the genetic dataset 
  # coords X and Y of each cell - number of cells - locus with 2 alleles
  # Args :
  #   rasK : raster of pop sizes
  #   dataSize : data that contains number of individuals in each cell
  #   nb_locus : number of locus wished
  #   initial_locus_value : genetic information for the first generation (parents)
  #   Option : sample/full_1col/2col_haploid/diploid
  #            sample : nind individuals are sampled from the full population
  #            full : number of individuals per cell equals carrying capacity
  #            1col/2col : for diploids, one columns or 2 columns representation
  #            haploid or diploid : haploid or diploid 
  #   nind = number of individuals sampled from the full carrying capacity population
  # Returns : a data frame that contains genetic data : x, y, Cell_numbers, 
  # 
  # Get coords for each cell
  coords = xyFromCell(object=rasK, cell=1:length(values(rasK[[1]])), spatial=FALSE)
  # Get a vector giving the number of the deme, repeted as many times as sampled individuals within (ex : 111144 )
  repet = switch(Option,
                 sample_1col_diploid = sort(rep(sample(rep(1:ncell(rasK),round(values(rasK))),nind,replace=TRUE),2)),
                 sample_2col_diploid = sample(rep(1:length(rasK),round(values(rasK))),nind,replace=TRUE),
                 sample_haploid = sample(rep(1:ncell(rasK),round(values(rasK))),nind,replace=TRUE),
                 full_1col_diploid = rep(1:length(rasK),round(values(rasK))*2),
                 full_2col_diploid = rep(1:length(rasK),round(values(rasK))),
                 full_haploid = rep(1:length(rasK),round(values(rasK))),            
  )
  # Make a data frame of genetic values :
  genes = switch(Option,
                 sample_1col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                 sample_2col_diploid  = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                 sample_haploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                 full_1col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                 full_2col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                 full_haploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),                 
  )
  # determine names of column, like "Locus 1","Locus 2" for 1_col_diploid and haploid, or "Locus_1.1", "Locus_1.2" for 2_col_diploid 
  colnames(genes) = switch(Option,
                           sample_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                           sample_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                           sample_haploid = paste("Locus",1:nb_locus, sep=""),
                           full_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                           full_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                           full_haploid = paste("Locus",1:nb_locus, sep=""),                 
  )
  # Merge sampling, geographic, genetic informations
  geneticData <- as.data.frame(coords[repet,])
  geneticData[,"Cell_numbers"] <- repet
  geneticData <- cbind(geneticData,genes) # add locus to geneticData 
  return(geneticData)
}