library(ape)
library(stringr)
library(markovchain)
library(matrixcalc)
library(MASS)

############## CLASS AND VALIDITY ####

setClass("TransitionForward",
         contains = "matrix",
         validity = function(object){
                        if (all(nrow(object)==0))stop("The matrix is empty.")
                        if (nrow(object)!=ncol(object))stop("The matrix is not square")
                      }
)

Demographic<-setClass("Demographic",
                      contains = "Landscape",
                      slots = c(K="numeric", R="numeric",TransiBackw="TransitionBackward",TransiBackwSym="TransitionBackward",TransiForw="TransitionForward"),
                      validity = function(object){
                        if(any(object@K<0))stop("K is negative")
                        if(any(object@R<0))stop("R is negative")
                      }
)

coalescentEvent <- setClass("coalescentEvent",
                          contains="environment",
                          slot = c(time="numeric",coalescing="integer",internal="integer",br_length="numeric",descendantTips="integer"),
                          prototype=list(time=10,coalescing=as.integer(c(1,2)),internal=as.integer(3),br_length=c(10,10),descendantTips=as.integer(c(1,2))),
                          )                          
                          
Coalescent <- setClass("Coalescent",
                       contains = "list",
                       validity = function(object){
                       if (any(unlist(lapply(object,function(x) class(x)!="coalescentEvent")))) 
stop("error in class Coalescent constructor : should be a list of coalescentEvent")
}
)


#Coalescent <- setClass("Coalescent",
#                       contains = "list",
#                       validity = function(object){
#                       if (!any(unlist(lapply(object,function(x) names(x)%in%c("time","coalescing","internal","br_length","descendantTips"))))) 
#stop("missing in coalescent list")
#}
#)


############## METHODS #####

setGeneric(
  name = "transitionMatrixBackward",
  def=function(object,model,symetric){return(standardGeneric("transitionMatrixBackward"))}
)

setGeneric(
  name = "transitionMatrixForward",
  def=function(param, meth){return(standardGeneric("transitionMatrixForward"))}
)


setGeneric(
  name = "createDemographic",
  def=function(object,model){return(standardGeneric("createDemographic"))}
)

setGeneric(
  name = "laplaceMatrix",
  def=function(object){return(standardGeneric("laplaceMatrix"))}
)

setGeneric(
  name = "ordinary_laplacian",
  def=function(object){return(standardGeneric("ordinary_laplacian"))}
)


setGeneric(
  name = "hitting_time_digraph",
  def=function(object){return(standardGeneric("hitting_time_digraph"))}
)

setGeneric(
  name = "commute_time_digraph",
  def=function(object){return(standardGeneric("commute_time_digraph"))}
)

setGeneric(
  name = "simul_coalescent",
  def=function(demographic,printCoal){return(standardGeneric("simul_coalescent"))}
)

setGeneric(
  name = "simul_multi_coal",
  def=function(demographic,printCoal,iteration){return(standardGeneric("simul_multi_coal"))}
)

setGeneric(
  name = "compare",
  def=function(demographic,popSize,printCoal,iteration){return(standardGeneric("compare"))}
)

setGeneric(
  name = "Collisionijk",
  def=function(Hitting_mat){return(standardGeneric("Collisionijk"))}
)

setGeneric(
  name = "linearizedFstDigraph",
  def=function(transition, popSize){return(standardGeneric("linearizedFstDigraph"))}
)

setGeneric(
  name = "coalescent_2_newick",
  def=function(coalescent){return(standardGeneric("coalescent_2_newick"))}
)

setGeneric(
  name = "linearizedFstUndigraph",
  def=function(transition, popSize){return(standardGeneric("linearizedFstUndigraph"))}
)

setGeneric(
  name = "commute_time_undigraph",
  def=function(object){return(standardGeneric("commute_time_undigraph"))}
)

setGeneric(
  name = "nodesInfo",
  def = function(coalescent){return(standardGeneric("nodesInfo"))}
)

setMethod(
  f ="[",
  signature = c(x="Demographic" ,i="character",j="missing"),
  definition = function (x ,i ,j , drop ){
    switch ( EXPR =i,
             "K" ={return(x@K)} ,
             "R" ={return(x@R)} ,
             "TransiBackw" ={return(x@TransiBackw)} ,
             "TransiBackwSym"={return(x@TransiBackwSym)},
             "TransiForw" = {return(x@TransiForw)},
             stop("This slots doesn't exist!")
    )
  }
)

setMethod(f="show",
          signature = "Demographic",
          definition = function(object){
            print(list(landscape=raster(object),K=object@K,R=object@R,TransiBackw=object@TransiBackw,TransiBackwSym=object@TransiBackwSym,TransiForw=object@TransiForw))
})

setMethod(f="transitionMatrixBackward",
          signature=c("Landscape","list","logical"),
          definition=function(object,model,symetric){
            if ((length(model$R)==1)&(length(model$K)==1)){transition = model$R * model$K * t(model$migration)}
            if ((length(model$R)>1)&(length(model$K)==1)){transition = t(matrix(model$R,nrow=length(model$R),ncol=length(model$R))) * model$K * t(model$migration)}
            if ((length(model$R)==1)&(length(model$K)>1)){transition = model$R * t(matrix(model$K,nrow=length(model$K),ncol=length(model$K))) * t(model$migration)}
            if ((length(model$R)>1)&(length(model$K)==1)){transition = t(matrix(model$R,nrow=length(model$R),ncol=length(model$R))) * lpar$K * t(model$migration)}
            if ((length(model$R)>1)&(length(model$K)>1)) {transition = t(matrix(model$R,nrow=length(model$R),ncol=length(model$R))) * t(matrix(model$K,nrow=length(model$K),ncol=length(model$K))) * t(model$migration)}
            if(symetric){
              symMat<-transition+t(transition)
              maxrow<-max(rowSums(symMat))
              diag(symMat)<-(sapply(1:nrow(symMat),function(i)maxrow-(rowSums(symMat)[i]-symMat[i,i])))
              t<-symMat/as.numeric(maxrow)
            }
            else{
              t<-transition/t(sapply(rowSums(transition),function(x)rep(x,ncol(transition))))
            }
            TransitionBackward(t)
          }
)

setMethod(
  f="transitionMatrixForward",
  signature=c("list","character"),
  definition=function(param, meth)
  {
    rs <- matrix(param$R,nrow=length(param$R),ncol=length(param$R))
    Ku = t(matrix(param$K,nrow=length(param$K),ncol=length(param$K)))
    leave = param$migration*(rs+1)*t(Ku); leave = leave - diag(leave)
    tMF<-switch (meth,
            non_overlap = param$migration*rs*Ku/colSums(rs * t(Ku) * param$migration),
            overlap = param$migration*(1+rs)*Ku/(colSums((1+rs)*t(Ku)*param$migration-t(leave))),
            stop("error in creation of transitionMatrixForward : the method does not exist !")
    )
    new(Class = "TransitionForward",tMF)
  }
)

setMethod(f="createDemographic",
          signature=c("Landscape","EnvDinModel"),
          definition=function(object,model){
            lpar<-runEnvDinModel(object,model)
            b<-transitionMatrixBackward(object,lpar,FALSE)
            bs<-transitionMatrixBackward(object,lpar,TRUE)
            f<-transitionMatrixForward(lpar,"non_overlap")
            new(Class = "Demographic",object,K=lpar$K,R=lpar$R,TransiBackw=b,TransiBackwSym=bs,TransiForw=f)
          }
)

setMethod(
  f = "nCellA",
  signature = "Demographic",
  definition = function(object){
    nCellA(object[[1]])
  }
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

setMethod(
  f="ordinary_laplacian",
  signature = "TransitionBackward",
  definition = function(object){
    markovB<-new("markovchain", states=dimnames(transition)[[1]], transitionMatrix=transition)
    PI<-diag(steadyStates(markovB)[1,])
    PI - PI%*%transition
  }
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

setMethod(
  f="commute_time_digraph",
  signature = "TransitionBackward",
  definition = function(object){
    mat<-hitting_time_digraph(object)
    sapply(1:ncol(mat),function(x)sapply(1:nrow(mat),function(y)mat[x,y]+mat[y,x]))
  }
)

          
          

############################################

setMethod("subset",
          signature="Coalescent",
          definition=function(x,internal=NULL,Time=NULL)
          {
            subCoal <- list()
            if (!is.null(internal)) 
            {
              if (!(internal%in%nodesInfo(x)$internals)) 
              {
                stop("subset is empty : internal node requested in subset does not exist in this genealogy") } else 
                {
                elements = grep(internal,lapply(x,function(x) slot(x,"internal"))) 
                subCoal[[1]] <- x[[internal]]
                tips <- nodesInfo(x)$tip
                i=1
                while (!all(subCoal[[i]]$coalescing%in%tips))
                  {
                  for (descendant in subCoal[[i]]$coalescing)
                    {
                    if (!(descendant%in%tips))
                      {
                      i=i+1
                      subCoal[[i]] <- x[[grep(descendant,lapply(x,function(x) x$internal))]]
                      }
                    }
                  }
                }
              subCoal <- lapply(order(unlist(lapply(subCoal,function(x) x$time))),function(i) subCoal[[i]])
              new("Coalescent",subCoal)
              } else {
              if (time>0) 
                {
                stop("time requested is negative in Coalescent subset") } else {
                Nodes <- which(unlist(lapply(x,function(x) x$time))<Time)
                tips <- nodesInfo(x)$tip
                subCoal
                i=0
                for (intNode in Nodes) 
                  {
                  i=i+1
                  subCoal[[i]] = x[[Nodes[intNode]]]
                  while (!all(subCoal[[i]]$coalescing%in%tips))
                    {
                    for (descendant in subCoal[[i]]$coalescing)
                      {
                      if (!(descendant%in%tips))
                        {
                        i=i+1
                        subCoal[[i]] <- x[[grep(descendant,lapply(x,function(x) x$internal))]]
                        }
                      }
                    }
                  }
                NewNodes <- lapply(subCoal,function(x) x$newNodes)
                uniqueElts <- !duplicated(NewNodes)
                subCoal <- lapply(uniqueElts,function(i) subCoal[[i]])
                subCoal <- lapply(order(unlist(lapply(subCoal,function(x) x$time))),function(i) subCoal[[i]])
                new("Coalescent",subCoal)
                }
             }
          }
)

setMethod(
  f="simul_coalescent",
  signature=c("Demographic","logical"),
  definition=function(demographic,printCoal)
  {
    prob_forward=NA
    N <- round(demographic["K"]);N[N==0]<-1
    coalescent = list()
    nodes = as.integer(rownames(xyFromCellA(demographic)));names(nodes)=as.character(nodes)
    cell_number_of_nodes <- as.integer(rownames(xyFromCellA(demographic)))             #point d'ou part la coalescent <- vecteurc numeric dont les valeur sont dans les cellules attribué
    names(cell_number_of_nodes) <- nodes
    parent_cell_number_of_nodes <- cell_number_of_nodes
    nodes_remaining_by_cell = list()
    time=0
    single_coalescence_events=0
    single_and_multiple_coalescence_events=0
    for (cell in 1:nCellA(demographic))
    {
      nodes_remaining_by_cell[[cell]] <- which(cell_number_of_nodes==cell)
    }
    while (length(unlist(nodes_remaining_by_cell))>1)
    {
      for (node in 1:length(parent_cell_number_of_nodes))
      {
        parent_cell_number_of_nodes[node] = sample(nCellA(demographic),size=1,prob=c(demographic["TransiBackw"][cell_number_of_nodes[node],]))
      }
      prob_forward[time] = sum(log(demographic["TransiForw"][parent_cell_number_of_nodes,cell_number_of_nodes]))
      time=time+1; if(printCoal==TRUE){if (round(time/10)*10==time) {print(time)}}
      for (cell in 1:nCellA(demographic))
      {
        nodes_remaining_by_cell[[cell]] <- as.integer(names(which(parent_cell_number_of_nodes==cell)))
        if (length(nodes_remaining_by_cell[[cell]])>1)
        {
          nbgenesremaining=length(nodes_remaining_by_cell[[cell]])
          smp = sample(N[cell],length(nodes_remaining_by_cell[[cell]]),replace=TRUE) # sample the parents among the N[cell] individuals of the cell
          parentoffspringmatrix <- matrix(smp,nrow=nbgenesremaining,ncol=N[cell])==matrix(1:N[cell],nrow=nbgenesremaining,ncol=N[cell],byrow=TRUE) # create a binary matrix of size offspring x parent that contains offspring parent relationships
          rownames(parentoffspringmatrix) <- nodes_remaining_by_cell[[cell]] # name the lines as the offspring nodes
          if (any(colSums(parentoffspringmatrix)>1) )
          {
            for (multiple in which(colSums(parentoffspringmatrix)>1))
            {
              single_coalescence_events = single_coalescence_events +1
              nodes_that_coalesce = names(which(parentoffspringmatrix[,multiple]))
              new_node <- max(nodes)+integer(1);nodes=nodes[!(names(nodes)%in%nodes_that_coalesce)];nodes=append(nodes,new_node);names(nodes)[length(nodes)]=new_node
              parent_cell_number_of_nodes <- append(parent_cell_number_of_nodes[!(names(parent_cell_number_of_nodes)%in%nodes_that_coalesce)],cell);names(parent_cell_number_of_nodes)[length(parent_cell_number_of_nodes)]<-new_node
              coalescent[[single_coalescence_events]] <- new("coalescentEvent",time=time,coalescing=as.integer(nodes_that_coalesce),internal=new_node,br_length=integer(),descendantTips=integer())
              nodes_remaining_by_cell[[cell]] <- append(nodes_remaining_by_cell[[cell]][!nodes_remaining_by_cell[[cell]] %in% nodes_that_coalesce],new_node)
              single_and_multiple_coalescence_events = single_and_multiple_coalescence_events + length(nodes_that_coalesce) - 1
            }
          }
        }
      }
      cell_number_of_nodes = parent_cell_number_of_nodes
    }
    coalescent <- new("Coalescent",coalescent)
    tips = NULL
    internals = NULL
    nodes = NULL
    times = NULL
    for (i in 1:length(coalescent))
    {
      nodes = append(nodes,c(coalescent[[i]]$coalescing,slot(coalescent[[i]],"internal")))
      internals = append(internals,coalescent[[i]]$internal)
      times = append(times,coalescent[[i]]$time)
    }
    nodes = as.numeric(levels(as.factor(c(nodes,internals))));nodes = nodes[order(nodes)]
    tips = nodes[!((nodes)%in%(internals))]
    # getting the branch length of each coalescing node
    # adding descendant tip information : internal nodes are attributed to descendant tip list
    for (i in 1:length(coalescent))#i=1
    {
      for (coalescing in coalescent[[i]]$coalescing)# coalescing = coalescent[[i]]$coalescing[1]
      {
        slot(coalescent[[i]],"tipDescendants") <- nodesInfo(subset(coalescent,internal=i))$tip
        if (coalescing %in% tips) {slot(coalescent[[i]],"br_length") <- append(slot(coalescent[[i]],"br_length"),coalescent[[i]]$time)
        } else {
          slot(coalescent[[i]],"br_length") <- slot(coalescent[[i]],"br_length")+slot(coalescent[[i]],"time")-times[which(internals==coalescing)]
        }
      }
    }
    list(coalescent=coalescent,prob_forward=sum(prob_forward))
  }
)


setMethod(
  f="simul_multi_coal",
  signature=c("Demographic","logical","numeric"),
  definition=function(demographic,printCoal,iteration){
    lapply(1:iteration,function(x) simul_coalescent(demographic,printCoal))
  }
)

setMethod("plot",
         signature="Coalescent",
         definition=function(x){
         plot.new()
         #tips=nodesInfo(object)$tip
         tips=gsub("\\(","",gsub("\\)","",gsub(";","",strsplit(coalescent_2_newick(object),",")[[1]])))
         coords = list()
         text(tips,x=((0:(length(tips)-1))/(length(tips)-1)),y=0)
         for (i in names(coalescent))
         {
           for (j in names(coalescent$coalescing))
           {
           coords[[j]] <- coalescent[[i]]
           }
         }
}
) 

setMethod(f="dist",
          signature="Coalescent",
          definition = function(x){
            ninfo = nodesInfo(x)
            coalPairMat <- diag(length(ninfo$nodes)) # we set the pairs of nodes that have coalesced together to the digaonal of tip nodes
            matDist = matrix(0,nrow=length(ninfo$tip),ncol=length(ninfo$tip), dimnames=list(ninfo$tip,ninfo$tip))
            # we set the matDist result to a 0 distance matrix
            coalesced=rep(FALSE,length(ninfo$tip))
            remaining=!coalesced
            remaining=ninfo$tip
            coalesced=NULL
            Time = 0 # we set to zero the time of the last coalescence
            # we loop over the coalescence events
            for (i in 1:length(x))
              {                
                # we add a distance to the pairs of nodes that have not coalesced that equals to the time since the last coalescence event 
                matDist = matDist+!coalPairMat*(x[[i]]$time - Time)  
                if length(coalPaiMat[x[[i]]$descendantTips==0) coalPaiMat[x[[i]]$descendantTips<-subset  
                coalPaiMat[x[[i]]$descendantTips,x[[i]]$descendantTips]<-TRUE
                matDist[remaining,ninfo$tip]=matDist[remaining,ninfo$tip]+slot(x[[i]],"br_length")
                coalesced = as.character(levels(as.factor(append(coalesced,remaining[remaining%in%slot(x[[i]],"coalescing")]))))
                remaining = remaining[!remaining%in%slot(x[[i]],"coalescing")]
              }
})

setMethod(
  f="compare",
  signature=c("Demographic","Landscape","logical","numeric"),
  definition=function(demographic,popSize,printCoal,iteration){
    coalescent<-simul_multi_coal(demographic,printCoal,iteration)
    lcoal<-lapply(1:iteration,function(n){
      dist(coalescent[[n]][[1]])
#      coal_2<-coalescent_2_newick(coalescent[[n]][[1]])
#      cat(coal_2, file = "ex.tre", sep = "\n")
#      tree<-read.tree("ex.tre")
#      cophenetic(tree)
    })
    a<-matrix(data = apply(sapply(lcoal,as.vector),1,mean),nrow = nrow(lcoal[[1]]),ncol = ncol(lcoal[[1]]))
    b<-linearizedFstDigraph(demographic["TransiBackw"],popSize)
    c<-linearizedFstUndigraph(demographic["TransiBackwSym"],popSize)
    d<-apply(popSize["distanceMatrix"],c(1,2),log)
    mat<-list(a,b,c,d)
    par(mfrow=c(2,2))
    for(i in 1:4){
      plot(bionj(mat[[i]]),main=title(switch(EXPR=as.character(i),
                                        "1"="Simul_coalescent",
                                        "2"="linearizedFstDigraph",
                                        "3"="linearizedFstUnDigraph",
                                        "4"="Stepping_Stone"))
           )
    }
    mat
  }
)


setMethod(
  f="Collisionijk",
  signature="matrix",
  definition=function(Hitting_mat)
  {
    Tijk=array(NA,dim=c(dim(Hitting_mat)[1],dim(Hitting_mat)[2],dim(Hitting_mat)[1]))
    for (k in 1:dim(Hitting_mat)[1]){
      for (i in 1:dim(Hitting_mat)[1]){
        for (j in 1:dim(Hitting_mat)[2]){
          Tijk[i,j,k] <- max(Hitting_mat[i,k],Hitting_mat[j,k])
        }
      }
    }
    Tijk
  }
)

############################################
setMethod(
  f = "valuesA",
  signature = "Landscape",
  definition = function(object){
    x=na.omit(values(object))
    colnames(x)=names(x)
    rownames(x) <- cellNumA(object)
    x
  }
)

setMethod(
  f="linearizedFstDigraph",
  signature=c("TransitionBackward","Landscape"),
  definition=function(transition, popSize)#popSize is raster class
  {
    H <- hitting_time_digraph(transition)
    dim2 <- dim(H);dim2[[3]]=2
    H2 <- array(c(H,t(H)),dim=dim2)
    MaxH <- apply(H2,c(1,2),max)
    genetic_dist = MaxH / (8*sum(valuesA(popSize))*nCellA(popSize))
    genetic_dist
  }
)

setMethod(
  f="coalescent_2_newick",
  signature="Coalescent",
  definition=function(coalescent)
  {
    tree=paste(" ",slot(coalescent[[length(coalescent)]],"internal")," ",sep="")
    for (i in length(coalescent):1)
    {
      Time = coalescent[[i]]$time
      coalesc <- as.character(slot(coalescent[[i]],"coalescing"))
      tree <- str_replace(tree,paste(" ",as.character(slot(coalescent[[i]],"internal"))," ",sep=""),paste(" ( ",paste(" ",coalesc," :",slot(coalescent[[i]],"br_length"),collapse=" ,",sep=""),") ",sep=""))
    }
    tree <- gsub(" ","",paste(tree,";",sep=""))
    tree
  }
)

setMethod(
  f="nodesInfo",
  signature="Coalescent",
  definition=function(coalescent)
  {
    nodes  = as.integer(levels(as.factor(unlist(lapply(coalescent, function(x) slot(x,"coalescing"))))))
    internals = unlist(lapply(coalescent, function(x) slot(x,"internal")))
    tip = nodes[!(nodes%in%internals)]
    x=list(tip=tip,internals=internals,nodes=nodes)
    x
  }
)

setMethod(
  f="linearizedFstUndigraph",
  signature=c("TransitionBackward","Landscape"),
  definition=function(transition, popSize)
  {
    commute_time <- commute_time_undigraph(transition)
    linearizedFst = commute_time / (16*sum(valuesA(popSize))*nCellA(popSize))
    linearizedFst
  }
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
    #mji = t(mij)
    commute_time = mii + mjj - 2*mij #- mji
    commute_time
  }
)
