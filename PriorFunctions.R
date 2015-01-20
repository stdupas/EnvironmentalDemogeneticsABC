uniform <- function(n, min, max){
  # A uniform function used to sample parameters in the specified prior distribution
  # 
  # Args:
  #   n: number of observations
  #   min: lower limit of the distribution
  #   max: upper limit of the distribution
  return(runif(n,min,max))
}