setGeneric(
  name = "laplaceMatrix",
  def=function(object){return(standardGeneric("laplaceMatrix"))}
)


setGeneric(
  name = "ordinary_laplacian",
  def=function(object){return(standardGeneric("ordinary_laplacian"))}
)

setGeneric(
  name = "commute_time_undigraph",
  def=function(object){return(standardGeneric("commute_time_undigraph"))}
)

setGeneric(
  name = "hitting_time_digraph",
  def=function(object){return(standardGeneric("hitting_time_digraph"))}
)

setGeneric(
  name = "genetDistUndigraph",
  def=function(object,popSize,mutation_rate){return(standardGeneric("genetDistUndigraph"))}
)

setGeneric(
  name = "genetDistDigraph",
  def=function(object,popSize,mutation_rate,method="Goldstein95"){return(standardGeneric("genetDistDigraph"))}
)

