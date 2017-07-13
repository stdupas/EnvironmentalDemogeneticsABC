library(ape)
library(stringr)
library(lattice)
library(markovchain)
library(matrixcalc)
library(abind)
library(rgdal)
library(raster)
library(MASS)
############## CREATION OF TransitionMatrix #######################################
r1<- raster(ncol=2, nrow=2)
r1[] <- rep(2:5,1)
r2<- raster(ncol=2, nrow=2)
r2[] <- rep(2,2:2)
s<- stack(x=c(r1,r2))
p1<-as.Date("2000-01-11")
variK<-c("l","t")
variR<-c("c","t")
paraK<-list(c(0,5),2)
paraR<-list(2,2)
reaK<-c(l="envelin",t="constant")
reaR<-c(l="constant",t="constant")
extent(s)<-c(0,2,0,2)
lscp1<-Landscape(rasterstack = s,period=p1,vars=vari)
modelK<-NicheModel(variables=vari,parameterList=para,reactNorms=rea)
modelR<-NicheModel(variables=vari,parameterList=para,reactNorms=rea)
m<-MigrationModel(shape="gaussian",param = (1/1.96))
edm1<-EnvDinModel(K=modelK,R=modelR,migration = m)
transi1<-createTransitionMatrix(lscp1,edm1)
xyFromCellA(lscp1)
###########

setGeneric(
  name = "laplaceMatrix",
  def=function(object){return(standardGeneric("laplaceMatrix"))}
)


setMethod(
  f = "laplaceMatrix",
  signature = "TransitionBackward",
  definition = function(object){
    matrixD = diag(rep(1,dim(object)[1])) # diagonal equals to 1
    laplacianMatrix = matrixD - object
    laplacianMatrix[is.na(laplacianMatrix)]<-0 # replace NA by 0
    #cat("laplacian",laplacianMatrix)
    return(laplacianMatrix)
  }
)

setGeneric(
  name = "ordinary_laplacian",
  def=function(object){return(standardGeneric("ordinary_laplacian"))}
)

setMethod(
  f="ordinary_laplacian",
  signature = "TransitionBackward",
  definition = function(object){
    markovB<-new("markovchain", states=dimnames(transition)[[1]], transitionMatrix=transition)
    PI<-diag(steadyStates(markovB)[1,])
    PI - PI%*%transition
  }
)


setGeneric(
  name = "commute_time_undigraph",
  def=function(object){return(standardGeneric("commute_time_undigraph"))}
)

setMethod(
  f="commute_time_undigraph",
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

setGeneric(
  name = "hitting_time_digraph",
  def=function(object){return(standardGeneric("hitting_time_digraph"))}
)

setMethod(
  f="hitting_time_digraph",
  signature = "TransitionBackward",
  definition = function(object){
    Ones <- rep(1,dim(object)[1])
    markovB<-new("markovchain", states=dimnames(object)[[1]], transitionMatrix=object)
    pi_<-steadyStates(markovB)[1,]
    PI <- diag(pi_)
    L <- PI - PI%*%object
    Z <- ginv(L + pi_%*%t(pi_))
    H <- Ones%*%t(diag(Z))-Z
    H
  }
)

setGeneric(
  name = "simul_coalescent",
  def=function(transitionList){return(standardGeneric("simul_coalescent"))}
)

setMethod(
  f="simul_coalescent",
  signature="list",
  definition=function(transitionList,geneticData)
  {
    prob_forward=NA
    N <- round(transitionList$K);N[N==0]<-1
    coalescent = list()
    nodes = as.numeric(rownames(geneticData));names(nodes)=as.character(nodes)
    cell_number_of_nodes <- geneticData[,"Cell_numbers"]
    names(cell_number_of_nodes) <- nodes
    parent_cell_number_of_nodes <- cell_number_of_nodes
    nodes_remaining_by_cell = list() 
    time=0 
    single_coalescence_events=0
    single_and_multiple_coalescence_events=0 
    for (cell in 1:ncellA(rasterStack))
    {
      nodes_remaining_by_cell[[cell]] <- which(cell_number_of_nodes==cell)
    }
    while (length(unlist(nodes_remaining_by_cell))>1) 
    {
      for (node in 1:length(parent_cell_number_of_nodes))
      {
        parent_cell_number_of_nodes[node] = sample(ncellA(rasterStack),size=1,prob=c(transitionList$backw[cell_number_of_nodes[node],]))
      }
      prob_forward[time] = sum(log(transitionList$forw[parent_cell_number_of_nodes,cell_number_of_nodes]))
      time=time+1; if (round(time/10)*10==time) {print(time)}
      for (cell in 1:ncellA(rasterStack))
      {
        nodes_remaining_in_the_cell = nodes_remaining_by_cell[[cell]] <- as.numeric(names(which(parent_cell_number_of_nodes==cell)))
      }
      prob_forward[time] = sum(log(transitionList$forw[parent_cell_number_of_nodes,cell_number_of_nodes]))
      time=time+1; if (round(time/10)*10==time) {print(time)}
      for (cell in 1:ncellA(rasterStack))
      {     
        nodes_remaining_in_the_cell = nodes_remaining_by_cell[[cell]] <- as.numeric(names(which(parent_cell_number_of_nodes==cell)))
        if (length(nodes_remaining_in_the_cell)>1) 
        {
          nbgenesremaining=length(nodes_remaining_in_the_cell)
          smp = sample(N[cell],length(nodes_remaining_in_the_cell),replace=TRUE)
          parentoffspringmatrix <- matrix(smp,nrow=nbgenesremaining,ncol=N[cell])==matrix(1:N[cell],nrow=nbgenesremaining,ncol=N[cell],byrow=TRUE)
          rownames(parentoffspringmatrix) <- nodes_remaining_in_the_cell
          if (any(colSums(parentoffspringmatrix)>1) )
          {
            for (multiple in which(colSums(parentoffspringmatrix)>1))
            {
              single_coalescence_events = single_coalescence_events +1
              nodes_that_coalesce = names(which(parentoffspringmatrix[,multiple]))
              new_node <- max(nodes)+1;nodes=nodes[!(names(nodes)%in%nodes_that_coalesce)];nodes=append(nodes,new_node);names(nodes)[length(nodes)]=new_node
              parent_cell_number_of_nodes <- append(parent_cell_number_of_nodes[!(names(parent_cell_number_of_nodes)%in%nodes_that_coalesce)],cell);names(parent_cell_number_of_nodes)[length(parent_cell_number_of_nodes)]<-new_node
              coalescent[[single_coalescence_events]] <- list(time=time,coalescing=as.numeric(nodes_that_coalesce),new_node=new_node)
              nodes_remaining_in_the_cell = nodes_remaining_by_cell[[cell]] <- append(nodes_remaining_in_the_cell[!nodes_remaining_in_the_cell %in% nodes_that_coalesce],new_node)
              single_and_multiple_coalescence_events = single_and_multiple_coalescence_events + length(nodes_that_coalesce) - 1
            }
          }
        }
      }
      cell_number_of_nodes = parent_cell_number_of_nodes
    }
    list(coalescent=coalescent,prob_forward=sum(prob_forward))
  }
)







############# manipulation #################
commute_time_undigraph(transi1)
hitting_time_digraph(transi1)

