transitionMatrixBackward <- function(r,K, migration){
  # Compute the backward probabilities of migration obtained with an isotropic migration hypothesis
  #
  # Args:
  #   r: a vector a values for growth rate for the various cells
  #   K: a vector of values for carrying capacities of the cells
  #   migration: a migration matrix
  #
  # Returns:
  #   The backward probabilities matrix
    
  if ((length(r)==1)&(length(K)==1)){transition = r * K * t(migration)}
  if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
  if ((length(r)==1)&(length(K)>1)){transition = r * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
  if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
  if ((length(r)>1)&(length(K)>1)) {transition = t(matrix(r,nrow=length(r),ncol=length(r))) * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
  
  transition = transition / rowSums(transition)  # Standardisation
  transition = transition * as.numeric(transition>1e-6) # removal of values below 1e-6
  transition = transition / rowSums(transition)  # Standardisation again
  transition
}


transitionMatrixForward <- function(r,K, migration, meth){
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