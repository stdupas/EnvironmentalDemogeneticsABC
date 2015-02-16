reproduction <- function(individuals, r)
{
  individuals$descendants <- rpois(length(individuals$cells),r[individuals$cells])
  individuals
}

migration <- function(M)
{
  
}
transitionDynamicsForward <- function(r,K, migration, meth){
  # Compute the forward probabilities of migration obtained with an isotropic hypothesis
  #
  # Args:
  #   r: a vector a values for growth rate for the various cells
  #   K: a vector of values for carrying capacities of the cells
  #   migration: a migration matrix
  #   meth: "non_overlap" or "overlap"
  #
  # Returns:
  #   The forward transition probabilities
  
  # note s the parental cell and u the descendant cell
  rs = matrix(r,nrow=length(r),ncol=length(r)) # growth rate in parental cell
  Ku = t(matrix(K,nrow=length(K),ncol=length(K))) # carrying capacity in descendant cell
  leave = migration*(1+rs)*t(Ku); leave = leave - diag(leave) # individuals that leave
  switch (meth,
          non_overlap = migration * rs * Ku / colSums(rs * t(Ku) * migration),
          overlap = migration * (1+rs) * Ku / (colSums((1+rs) * t(Ku) * migration - t(leave)))
  )
}

repnDispMutFunction <- function(tipDemes,transitionForward, K){
  # Calcul for reproduction and dispersion 
  # random choice of individuals and at the same time of their target cells in the transition matrix
  
  ncell_transition <- dimGeneticData[1]*dim(transitionmatrice)[1]
  transition_celnu <- sample(ncell_transition,dimGeneticData[1],replace=FALSE,transitionmatrice[geneticData[,"Cell_numbers"],])
  transition_col <- ceiling(transition_celnu/dimGeneticData[1])
  transition_line <- transition_celnu%%(dimGeneticData[1]);transition_line[transition_line==0]<-dimGeneticData[1]
  cell_numbers_sampled <- geneticData[transition_line,"Cell_numbers"]
  geneticData <- geneticData[transition_line,]
  geneticData[,"Cell_numbers"] <- transition_col
  locusCols = grep("Locus", colnames(geneticData))
  step = 2
  mu = mutationRate # mutation rate
  liability = runif(prod(dimGeneticData), 0, 1) # mutation liability
  liability = as.data.frame(matrix(liability, ncol = length(locusCols), nrow = dimGeneticData[1]))
  geneticData[,locusCols] = geneticData[,locusCols] + ((liability<mu/2)*step - (liability>(1-mu/2))*step) 
  #print(c("mutrat",(sum(liability<mu/2)+sum(liability>(1-mu/2)))/(length(grep("Locus", colnames(geneticData)))*dimGeneticData[1])))
  geneticData
}
