

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
  f="genetDistUndigraph",
  signature = "TransitionBackward",
  definition = function(object,popSize,mutation_rate){
    commute_time <- commute_time_undigraph(object)
    #genetic_dist = commute_time / (8* popSize)
    #genetic_dist = commute_time / (8* (sum(popSize)/(dim(popSize)[1]*dim(popSize)[2])))
    genetDist =  (commute_time/4+2*sum(valuesA(popSize))) * mutation_rate
    genetDist
  }
)


setMethod(
  f="genetDistDigraph",
  signature = "TransitionBackward",
  definition = function(object,popSize,mutation_rate,method="Goldstein95"){
    H <- hitting_time_digraph(object)
    dim2 <- dim(H);dim2[[3]]=2
    H2 <- array(c(H,t(H)),dim=dim2)
    MinH <- apply(H2,c(1,2),min)
    #H2[,,1] <- MinH; H2[,,2] = (H+t(H))/2
    #MinH2 <- apply(H2,c(1,2),min)
    genetic_dist = (MinH/2 + 2*sum(valuesA(popSize)))* mutation_rate
    genetic_dist
  }
)


setMethod(
  f="simul_coal_genet",
  signature="",
  definition=function(geneticData,rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp,mutation_rate=1E-2,initial_locus_value=200,mutation_model="stepwise",stepvalue=2,locusnames=NA)
  {
    if (is.na(locusnames)) {
      locusnames <- grep("ocus",colnames(geneticData),value=TRUE)
      if (sum(grep("\\.1",locusnames,value=FALSE)|
              grep("\\.2",locusnames,value=FALSE))==length(locusnames)) {
        allelenames<-locusnames} else {
          allelenames <- paste(rep(locusnames,each=2),".",1:2,sep="")
        }                           
    }
    tip_genotype <- geneticData[,1:3]
    coalescent <- list()
    transitionMatrixList <- transitionMatrice(rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp)
    for (allele in allelenames){
      coalescent[[allele]] <- simul_coalescent(transitionMatrixList,geneticData)
      coalescent_mut=add_br_length_and_mutation(coalescent[[allele]]$coalescent,mutation_rate,initial_locus_value,allele)
      coaltable <- coalist_2_coaltable(coalescent_mut,allele,initial_locus_value)
      genetic_values=genetics_of_coaltable(coaltable,initial_locus_value,mutation_model,stepvalue)
      rownames(genetic_values) <- genetic_values$coalescing
      tip_genotype[,allele]<-genetic_values[as.character(1:dim(geneticData)[1]),"genotype"]
    }
    list(coalescent=coalescent,mutation_rate=mutation_rate,tip_genotype=tip_genotype,
         transitionAndK=transitionMatrice(rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp))
  }
)


setMethod(
  f="simul_coalescent",
  signature="",
  definition=function(transitionList,geneticData)
  {
    prob_forward=NA
    N <- round(transitionList$K);N[N==0]<-1
    coalescent = list()
    nodes = as.numeric(rownames(geneticData));names(nodes)=as.character(nodes)
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


setMethod(
  f="CreateGenetArray",
  signature="",
  definition=function(rasK, nb_locus, initial_locus_value,Option="sample_1_col_diploid",nind="Ne")
  {
    if (nind%in%c("Ne","Ne2","Ne10")) nind <- switch(nind,
                                                     Ne = floor(sum(valuesA(rasK))),
                                                     Ne2= 2*floor(sum(valuesA(rasK))),
                                                     Ne10= 10*floor(sum(valuesA(rasK)))
    )
    coords = xyFromCellA(rasK)
    repet = switch(Option,
                   sample_1col_diploid= sort(rep(sample(rep(1:ncellA(rasK),round(valuesA(rasK))),nind,replace=TRUE),2)), 
                   sample_2col_diploid= sample(rep(1:length(rasK),round(valuesA(rasK))),nind,replace=TRUE), 
                   sample_haploid= sample(rep(1:ncellA(rasK),round(valuesA(rasK))),nind,replace=TRUE),
                   full_1col_diploid= rep(1:length(rasK),round(valuesA(rasK))*2), 
                   full_2col_diploid= rep(1:length(rasK),round(valuesA(rasK))), 
                   full_haploid= rep(1:length(rasK),round(valuesA(rasK))),
                   homogenous_1col= rep(1:length(rasK),each=nind/length(rasK)*2),
                   homoNonEmtpyCells1col= rep((1:length(rasK))[which(round(valuesA(rasK))>1)],each=nind/length(rasK)*2)
    )
    geneticData <- as.data.frame(coords[repet,]) ; geneticData[,"Cell_numbers"] <- repet
    genes = switch(Option,
                   sample_1col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                   sample_2col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                   sample_haploid =      as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                   full_1col_diploid =   as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                   full_2col_diploid =   as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                   full_haploid =        as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                   homogenous_1col =     as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                   homogenous_2col =     as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus)*2)),
                   homoNonEmtpyCells1col=as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus)))
    )
    colnames(genes) = switch(Option,
                             sample_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                             sample_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                             sample_haploid = paste("Locus",1:nb_locus, sep=""),
                             full_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                             full_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                             full_haploid = paste("Locus",1:nb_locus, sep=""),                 
                             homogenous_1col = paste("Locus",1:nb_locus, sep=""),                 
                             homoNonEmtpyCells1col= paste("Locus",1:nb_locus, sep=""),                 
                             homogenous_2col =  paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep="")
    )
    geneticData
  }
)


############################ METHODS NON ADAPTÉE #############################################

setMethod(
  f="ssr",
  signature="",
  definition=function(p)
  {
    popSize = K_Function(rasterStack = donneesEnvironmentObs, p[1], p[2:(nblayers+1)])
    matrice_migration = migrationMatrix(donneesEnvironmentObs,shapeDisp, p[(nblayers+2):(nbpDisp+nblayers+1)])
    matrice_transition = transitionMatrixBackward(popSize, matrice_migration)
    matrice_laplace = laplaceMatrix(matrice_transition)
    commute_time = resistDist(matrice_laplace)
    linearizedFst_att = linearizedFstUndigraph(commute_time, p[(nbpDisp+nblayers+2)])
    linearizedFst_att = linearizedFst_att[finalGenetData$Cell_Number,finalGenetData$Cell_Number]
    return(mean((linearizedFst_obs - linearizedFst_att)^2))
  }
  
  
)


setMethod(
  f="test_stabilite_a_value",
  signature="",
  definition=function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice)
  {
    vecteur_a_value <-c(0)
    GenDist=list()
    print(i)
    geneticData <- repnDispMutFunction(geneticData, dimGeneticData, mutationRate, transitionmatrice) 
    Genotypes = geneticData[,grep("Locus", colnames(geneticData), fixed = T)]
    popmbrship=geneticData[,"Cell_numbers"]
    Fst = Fstat(Genotypes,nCell,popmbrship,ploidy=2)
    Fis = Fst[[1]];Fst=Fst[[2]]
    GenDist[[i]] <- Fst/(1-Fst)
    if((i>90) && (i%%30 == 0)){
      if(var(GenDist[[(i-30):i]])> var(GenDist[[(i-60):(i-30)]]) 
         && var(GenDist[[(i-30):i]])> var(GenDist[[(i-90):(i-60)]])){
        return(list(geneticData, GenDist))
        break
      }
    }
  }
)


setMethod(
  f="get_ultrametric_nj_tree",
  signature="",
  definition=function(tip_genotype,mutation_model,step_value)
  {
    nj_tree <- get_nj_tree(tip_genotype,mutation_model,step_value)
    chronos(nj_tree)
  }
)


setMethod(
  f="plot_coalescent",
  signature="",
  definition=function(coalescent_simulated,with_landscape=FALSE,rasK=NULL,legend_right_move=-.2,file=NA)
  {
    coalescent <- coalescent_simulated$coalescent
    gDat <- coalescent_simulated$tip_genotypes
    if (!is.na(file)) pdf(file)
    if (with_landscape) {par(mfrow=c(3,1),oma=c(0,0,0,4),xpd=TRUE)}else{par(mfrow=c(2,1),oma=c(0,0,0,4),xpd=TRUE)}
    tipcells <- gDat$Cell_numbers[as.numeric(coalescent_2_phylog(coalescent)$tip.label)]
    tipcols = rainbow(ncellA(rasK))[tipcells]
    plot(coalescent_2_phylog(coalescent),direction="downward",tip.color=tipcols)
    legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(legend_right_move,0))
    text(x=7,y=0,paste("tmra = ",tmra(coalescent)))
    njtree<-get_nj_tree(coalescent,mutation_model)
    plot(njtree,direction="downward")
    if (with_landscape) {plot(rasK)}
    if (!is.na(file)) dev.off()
  }
  #
  #
  
)


setMethod(
  f="coalescenceProb",
  signature="",
  definition=function(Qlist,time)
  {
    Ndeme <- (2*dim(Qlist[[1]])[1]+1/4)^.5-.5
    cTPD <- array(NA,dim=c(Ndeme,Ndeme,1))
    condition=TRUE
    t=1
    QPower <- Q <- Qlist[["Q"]]
    cTPD[,,t] <- matrix(colSums((matrix(1,nrow(QPower),ncol(QPower))-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme)
    while (condition){
      t=t+1
      QPowerBefore <- QPower
      QPower <- Q%*%QPowerBefore
      cTPD <- abind(cTPD,matrix(colSums((QPowerBefore-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme),along=3)
      condition=!all(colSums(cTPD,3)>0.9999)
    }
    cTPD
  }
  
)


setMethod(
  f="CreateGenetArray",
  signature="",
  definition=function(rasK, nb_locus, initial_locus_value,Option="sample_1_col_diploid",nind="Ne")
  {
    if (nind%in%c("Ne","Ne2","Ne10")) nind <- switch(nind,
                                                     Ne = floor(sum(valuesA(rasK))),
                                                     Ne2= 2*floor(sum(valuesA(rasK))),
                                                     Ne10= 10*floor(sum(valuesA(rasK)))
    )
    coords = xyFromCellA(rasK)
    repet = switch(Option,
                   sample_1col_diploid= sort(rep(sample(rep(1:ncellA(rasK),round(valuesA(rasK))),nind,replace=TRUE),2)), 
                   sample_2col_diploid= sample(rep(1:length(rasK),round(valuesA(rasK))),nind,replace=TRUE), 
                   sample_haploid= sample(rep(1:ncellA(rasK),round(valuesA(rasK))),nind,replace=TRUE),
                   full_1col_diploid= rep(1:length(rasK),round(valuesA(rasK))*2), 
                   full_2col_diploid= rep(1:length(rasK),round(valuesA(rasK))), 
                   full_haploid= rep(1:length(rasK),round(valuesA(rasK))),
                   homogenous_1col= rep(1:length(rasK),each=nind/length(rasK)*2),
                   homoNonEmtpyCells1col= rep((1:length(rasK))[which(round(valuesA(rasK))>1)],each=nind/length(rasK)*2)
    )
    geneticData <- as.data.frame(coords[repet,]) ; geneticData[,"Cell_numbers"] <- repet
    genes = switch(Option,
                   sample_1col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                   sample_2col_diploid = as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                   sample_haploid =      as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                   full_1col_diploid =   as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))), 
                   full_2col_diploid =   as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus*2))), 
                   full_haploid =        as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                   homogenous_1col =     as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus))),
                   homogenous_2col =     as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus)*2)),
                   homoNonEmtpyCells1col=as.data.frame(array(initial_locus_value,dim=c(length(repet),nb_locus)))
    )
    colnames(genes) = switch(Option,
                             sample_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                             sample_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                             sample_haploid = paste("Locus",1:nb_locus, sep=""),
                             full_1col_diploid = paste("Locus",1:nb_locus, sep=""), 
                             full_2col_diploid = paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep=""), 
                             full_haploid = paste("Locus",1:nb_locus, sep=""),                 
                             homogenous_1col = paste("Locus",1:nb_locus, sep=""),                 
                             homoNonEmtpyCells1col= paste("Locus",1:nb_locus, sep=""),                 
                             homogenous_2col =  paste("Locus",sort(rep(1:nb_locus,2)),".", 1:2, sep="")
    )
    geneticData
  }
)
  
  
setMethod(
  f="plotGeneticData",
  signature = "",
  definition = function(geneticData, EnvironmentalDataObserved)
  {
    colnames(geneticData)[1:2] <- c("x","y")
    geneticData = SpatialPixelsDataFrame(points = geneticData[,c("x","y")], data = geneticData[,])
    plot(EnvironmentalDataObserved[[1]])
    plot(geneticData, add = T)
  }
)


setMethod(
  f="linearizedFstDigraph",
  signature="",
  definition=function(transition, popSize)#popSize is raster class
  {
    H <- hitting_time_digraph(transition)
    dim2 <- dim(H);dim2[[3]]=2
    H2 <- array(c(H,t(H)),dim=dim2)
    MaxH <- apply(H2,c(1,2),max)
    genetic_dist = MaxH / (8*sum(valuesA(popSize))*ncellA(popSize))
    genetic_dist
  }
  
)


setMethod(
  f="genetDist",
  signature="",
  definition=function(tip_genotype,method="deltaMu",stepvalue=2,byCell=TRUE)
  {
    genetDist =switch(method,
                      deltaMu=deltaMu(tip_genotype,stepvalue=stepvalue))
    genetDist
  }
  
)


setMethod(
  f="OneCol23Dims",
  signature="",
  definition=function(gDat)
  {
    first <- gDat[2*(1:(dim(gDat)[1]/2))-1,grep("ocus",colnames(gDat))]  
    second <- gDat[2*(1:(dim(gDat)[1]/2)),grep("ocus",colnames(gDat))]  
    dimN <- dimnames(first);dimN[[3]]<-c("al.1","al.2")
    array(c(unlist(first),unlist(second)),dim=c(dim(first),2),dimnames=dimN)
  }
  
)


setMethod(
  f="summary_stat",
  signature="",
  definition=function(geneticDataObs,geneticDataSimulList,log_lik_simul_list)
  {
    
  }
  
  
  
  
  #
  #
  #
  
  
)


setMethod(
  f="SSw",
  signature="",
  definition=function(gDat)
  {
    2*(1-Qwithin_pop(gDat)) 
  }
  
)


setMethod(
  f="",
  signature="",
  definition=#
    
    
    
    
    
    
    
    
    
)


setMethod(
  f="populationSize",
  signature="",
  definition=function(donneesEnvironmentObs, p, shapes)
  {
    populationSize <- donneesEnvironmentObs[[1]]
    populationSize[cellNumA(populationSize)] <- ReactNorm(valuesA(donneesEnvironmentObs), p, shapes)[1,]
    populationSize
  }
  
)


setMethod(
  f="R_Function",
  signature="",
  definition=function(rasterStack, alpha, beta)
  {
    if(nlayers(rasterStack)>1){
    }
    else{ R = exp(as.matrix(alpha+beta*rasterStack)) }
    R
  }
  
  {
    pairs = NULL;i=0
    while (i < dim(p)[2])
    {
      i=i+1;j=i+1
      first= colnames(p)[i]
      while (j <= dim(p)[2])
      {
        second = colnames(p)[j]
        j=j+1
        pairs=rbind(pairs,c(first,second))
      }
    }
    pairs
    dev.off()
    par(mfrow=c(1,1))
    Data= as.data.frame(Data)
    for (i in 1:dim(pairs)[1])
    {
      print(i)
      form = as.formula(paste("z~",paste(pairs[i,],collapse="*"),sep=""))
      Data[,"z"]=ReactNorm(Data,p,shapes)[,"Y"]
    }
  }
  
  #
  #
)


setMethod(
  f="envelinear",
  signature="",
  definition=function(X,p,log=FALSE)
  {
    Yxmin = p[rep("Yxmax",dim(X)[1]),colnames(X)]
    Yxmax = p[rep("Yxmin",dim(X)[1]),colnames(X)]
    Xmin=  p[rep("Xmin",dim(X)[1]),colnames(X)]
    Xmax=  p[rep("Xmax",dim(X)[1]),colnames(X)]
    a = (Yxmin - Yxmax) / (Xmin - Xmax)
    b = Yxmin - Xmin * a
    if (log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
  }
  
)


setMethod(
  f="tmra",
  signature="",
  definition=function(coalescent_simulated)
  {
    coalescent_simulated$coalescent[[length(coalescent_simulated$coalescent)]]$time
  }
  #
  
)


setMethod(
  f="linear",
  signature="",
  definition=function(X,p,Log=FALSE)
  {
    Yx1 = p[rep("Yx1",dim(X)[1]),colnames(X)]
    Yx0 = p[rep("Yx0",dim(X)[1]),colnames(X)]
    a = (Yx1 - Yx0)
    b = Yx0
    if (Log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
  }
  
  
  
  
)


setMethod(
  f="probgenet",
  signature="",
  definition=function(transition,gDat)
  {
    if (all(unlist(strsplit(colnames(gDat)[grep("ocus",colnames(gDat))],"\\."))
            [1:length(colnames(gDat)[grep("ocus",colnames(gDat))])*2]==
            rep(c("1","2"),length(colnames(gDat)[grep("ocus",colnames(gDat))])))){
            gDat <- TwoCols2OneCol(gDat)
  } 
  for (allele in grep("ocus",colnames(gData))){
    
  }
  }




    )
  
  
  setMethod(
    f="Aggregate_and_adjust_raster_to_data",
    signature="",
    definition=function(Envir_raster_stack,release,recovery,extend_band_size,aggregate_index)
    {
      samples <- SpatialPoints(rbind(na.omit(release[,c("X","Y")]),na.omit(recovery[,c("X","Y")])))
      if (aggregate_index > 1) { Envir_raster_stack <- aggregate(crop(Envir_raster_stack,extent(samples)+extend_band_size), fact=aggregate_index, fun=mean, expand=TRUE, na.rm=TRUE) } else {
        Envir_raster_stack <- crop(Envir_raster_stack,extent(samples)+extend_band_size)
      }
      Envir_raster_stack
    }
    
    #
    
    
    
  )
  
  
  setMethod(
    f="mmute",
    signature="",
    definition=function(cells=c(1,2),transitionmatrice)
    {
      tmpcell <- cells[1];t=1
      while (tmpcell != cells[2])
      {
        tmpcell =  sample(dim(transitionmatrice)[2],size=1,prob=c(transitionmatrice[tmpcell,]))
        t=t+1
      }
      hit=TRUE
      while (tmpcell != cells[1])
      {
        tmpcell =  sample(dim(transitionmatrice)[2],size=1,prob=c(transitionmatrice[tmpcell,]))
        t=t+1
      }
      commute=TRUE
      t
    }
  )
  
  
  setMethod(
    f="coalescence_prob_time_distribution_matrix",
    signature="",
    definition=function(transition,max_time_interval=4,rasK,threshold=1E-6)
    {
      popSizes_receiving <- matrix(valuesA(rasK),nrow=dim(transition)[1],ncol=dim(transition)[2],byrow=TRUE)
      coal_column <- 1/(2*popSizes_receiving)
      coal_column[coal_column>1]<-1
      coal_column_tiMax <- 1-(1-coal_column)^max_time_interval
      not_coalesced_prob <- array(NA, dim=c(dim(transition),1))
      time_points <- 1
      occurence_prob[,,1] <-  diag(1,dim(transition)[1])
      occurence_and_coalescence[,,1] <- occurence_prob[,,1] / (2*popSizes_receiving)
      coalescence_prob[,,1] <- occurence_prob[,,1] %*% t(occurence_and_coalescence[,,1])
      coalescence_prob[,,1][coalescence_prob[,,1]>1]=1
      not_coalesced_prob[,,1] <- 1-coalescence_prob[,,1]
      condition=TRUE
      i=1
      maxtransition <- matrix.power(transition,max_time_interval)
      while (condition){
        i=i+1
        time_interval <- time_points[i-1]
        if (time_interval>max_time_interval) {
          time_interval <- max_time_interval;transition_powered = maxtransition
          coal_column_ti <- coal_column_tiMax
        } else {
          transition_powered = matrix.power(transition,time_interval)
          coal_column_ti <- 1-(1-coal_column)^time_interval
        }
        time_points[i] <- time_points[i-1]+time_interval
        occurence_prob <- abind(occurence_prob,occurence_prob[,,i-1] %*% transition_powered,along=3)
        occurence_and_coalescence <- abind(occurence_and_coalescence,occurence_prob[,,i] *coal_column_ti,along=3)
        coalescence_prob <- abind(coalescence_prob,occurence_prob[,,i] %*% t(occurence_and_coalescence[,,i]),along=3)
        coalescence_prob[,,i][coalescence_prob[,,i]>1]=1
        not_coalesced_prob <- abind(not_coalesced_prob,not_coalesced_prob[,,i-1]*(1-coalescence_prob[,,i]),along=3)
        coalescence_prob[,,i] <- coalescence_prob[,,i]*not_coalesced_prob[,,i-1]
        condition = !all(not_coalesced_prob[,,i]<threshold)
      }
      dimnames(coalescence_prob) <- list(1:dim(transition)[1],1:dim(transition)[2],time_points)
      expected_coalescence_times <- matrix(NA,nrow=dim(transition)[1],ncol=dim(transition)[2])
      for (i in 1:dim(transition)[1]){
        for (j in 1:dim(transition)[2]){
          expected_coalescence_times[i,j] <- coalescence_prob[i,j,]%*%time_points
        }
      }
      list(coalescent_prob=coalescence_prob,exp_times=expected_coalescence_times)
    }
    
    
    
    
    
    
    #
    
  )
  
  
  setMethod(
    f="hitting_time_digraph",
    signature="",
    definition=function(transition)
    {
      Ones <- rep(1,dim(transition)[1])
      markovB<-new("markovchain", states=dimnames(transition)[[1]], transitionMatrix=transition)
      pi_<-steadyStates(markovB)[1,]
      PI <- diag(pi_)
      L <- PI - PI%*%transition
      Z <- ginv(L + pi_%*%t(pi_))
      H <- Ones%*%t(diag(Z))-Z
      H
    }
    
    
  )
  
  
  setMethod(
    f="id_mat_prob",
    signature="",
    definition=function(A,B)
    {
      id_mat_prob <- array(NA,dim=c(dim(A)[1],dim(A)[1]),dimnames=list(unlist(dimnames(A)[1]),unlist(dimnames(A)[1])))
      for (i in 1:dim(A)[1]){
        for (j in 1:dim(A)[1]){
          id_mat_prob[i,j] <- mean(A[i,]==B[j,])
        }
      }
      id_mat_prob
    }
    
  )
  
  
  setMethod(
    f="genetics_of_coaltable",
    signature="",
    definition=function(coaltable,initial_locus_value,mutation_model="stepwise",stepvalue=2)
    {
      switch(mutation_model,
             stepwise = stepwise(coaltable,initial_locus_value,stepvalue,locusnames),
             my_ass = "rien"
      ) 
      #stepwise(coaltable,initial_locus_value,stepvalue, locusnames)
    }
    
    #
    #
  )
  
  
  setMethod(
    f="gridRepnDispFunction",
    signature="",
    definition=function(dynamics,r,K,d=.9,ptG, migration,overlapping=TRUE)
    {
      Nt = values(dynamics)[,dim(values(dynamics))[2]]*(1-d) + r*N*(K-N/K)
      esperance[K==0] <- 0  
    }
    
  )
  
  
  setMethod(
    f="K_Function",
    signature="",
    definition=function(rasterStack, p, shapes)
    {
      ReactNorm(valuesA(rasterStack),p,shapes)
    }
    
  )
  
  
  setMethod(
    f="coalescent_2_phylog",
    signature="",
    definition=function(coalescent)
    {
      read.tree(text=coalescent_2_newick(coalescent))
    }
    
    #
    #
    
  )
  
  
  setMethod(
    f="deme_coocurence_probability",
    signature="",
    definition=function(pGenes,transition,time)
    { 
      p=pGene1%*%t(pGene2)
    }
    
  )
  
  
  setMethod(
    f="forward_simul_landpopsize",
    signature="",
    definition=function(N0,p, migration)
    {
      
    }
    
  )
  
  
  setMethod(
    f="linearizedFstDigraph1",
    signature="",
    definition=function(transition, popSize)#popSize is raster class
    {
      H <- hitting_time_digraph(transition)
      Coalij <- H
      Coal = array (NA,dim=list(dim(H)[1],dim(H)[2],dim(H)[2]))
      for (i in 1:dim(Coal)[1]){
        for (j in 1:dim(Coal)[2]){
          for (k in 1:dim(Coal)[3]){
            Coal[i,j,k] <- max(H[i,k],H[j,k])
          }
          Coalij[i,j] <- mean(Coal[i,j,])
        }
      }
      linearizedFst = Coalij / (2*sum(valuesA(popSize))*ncellA(popSize))
      linearizedFst
    }
    
  )
  
  
  setMethod(
    f="comuteTimeDigraph",
    signature="",
    definition=function(transition, popSize)#popSize is raster class
    {
      H <- hitting_time_digraph(transition)
      commuteTime = (H+t(H)) 
    }
    
    #
    #
    #
    #
    
  )
  
  
 
  
  
  setMethod(
    f="S",
    signature="",
    definition=function(rasterStack)
    {
      S=matrix(0,nrow=ncellA(rasterStack),ncol=ncellA(rasterStack))
    d=dim(rasterStack)[1:2]
    ij <- t(t(xyFromCellA(rasterStack))/(res(rasterStack))+.5)
    for (donor in 1:ncellA(rasterStack)){
      for (receiver in 1:ncellA(rasterStack)){
        i <- abs(ij[donor,1]-ij[receiver,1])
        j <- abs(ij[donor,2]-ij[receiver,2])
        for (k in 0:(d[1]-1)){
          for (l in 0:(d[2]-1)){
            S[donor,receiver]=S[donor,receiver]+(2-as.numeric(k==0))*(2-as.numeric(l==0))*(1-cos(pi*i*k/d[1])*cos(pi*j*l/d[2]))/(1-.5*(cos(pi*k/d[1])+cos(pi*l/d[2])))
          }
        }
      }
    }
    S
    }
  )
  
  
  setMethod(
    f="OneCol2TwoCols",
    signature="",
    definition=function(gDat)
    {
      NonGeneticCols <- gDat[,-grep("ocus",colnames(gDat))]
      gDat <- gDat[,grep("ocus",colnames(gDat))]
      if (dim(gDat)[1] %% 2 != 0) stop("odd number of rows in suposedly one col genetic table")
      pair <- (1:(dim(gDat)[1]/2))*2;impair=pair-1
      first <- gDat[impair,]
      colnames(first) <- paste(colnames(first),"1",sep=".")
      second <- gDat[pair,]
      colnames(second) <- paste(colnames(second),"2",sep=".")
      Twocols<-cbind(first,second)
      neworder <- rep(c(0,dim(gDat)[2]),dim(gDat)[2])+rep(1:dim(gDat)[2],each=2)
      Twocols <- Twocols[,neworder]
      cbind(NonGeneticCols[impair,],Twocols)
    }
    
  )
  
  
  setMethod(
    f="ordinary_laplacian",
    signature="",
    definition=function(transition)
    {
      markovB<-new("markovchain", states=dimnames(transition)[[1]], transitionMatrix=transition)
      PI<-diag(steadyStates(markovB)[1,])
      PI - PI%*%transition
    }
    
    
    
  )
  
  
  setMethod(
    f="likelihoodCoalescent",
    signature="",
    definition=function(coalescent)
    {
      
    }
    
  )
  
  
  setMethod(
    f="simulationGenet",
    signature="",
    definition=function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, mutationRate, nbLocus, initial_locus_value, shapeDisp, pDisp, nb_generations,ind_per_cell=30)
    {
      Rast_K <- donneesEnvironmentObs ; valuesA(Rast_K) <- K
      Rast_r <- donneesEnvironmentObs ; valuesA(Rast_r) <- r
      geneticData = CreateGenetArray(Rast_K, nbLocus,initial_locus_value,Option="full_pop")
      geneticData[,"Cell_number_init"] <- geneticData[,"Cell_numbers"]
      dimGeneticData = dim(geneticData)
      migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
      transitionmatrice = transitionMatrixBackward(Npop = K, migration= migrationM)
      geneticDataFinal = test_stabilite_a_value(geneticData, mutationRate, dimGeneticData, nb_generations,transitionmatrice)
      a_value_obs = geneticDataFinal[[2]]
      geneticDataFinal = geneticDataFinal[[1]]
      return(list(geneticDataFinal, a_value_obs))
    }
    
    #
    #
    
  )
  
  
  setMethod(
    f="get_nj_tree",
    signature="",
    definition=function(tip_genotype,mutation_model,step_value)
    {
      if (mutation_model=="stepwise") distmat <- deltaMu(tip_genotype,stepvalue=step_value)
      nj(as.matrix(distmat))
    }
    
  )
  
  
  setMethod(
    f="laplaceMatrix",
    signature="",
    definition=function(transitionMatrix)
    {
      laplacianMatrix = matrixD - transitionMatrix
      laplacianMatrix
    }
    #
    
  )
  
  
  setMethod(
    f="genetDistAbsorbingMethod",
    signature="",
    definition=function(transition,N,mutation_rate)
    {
      Qlist<-absorbingTransition(transition,N)
      absorbingTime(fundamentalMatrixAbsorbingChain(Qlist$Q),Qlist$Qline)*mutation_rate
    }
    
    #
    #
    #
    
  )
  
  
  setMethod(
    f="proportional",
    signature="",
    definition=function(X,p,Log=FALSE)
    {
      if (Log) {log(p[rep("a",dim(X)[1]),colnames(X)]*X)} else {
        p[rep("a",dim(X)[1]),colnames(X)]*X
      }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  )
  
  
  setMethod(
    f="linearizedFstUndigraph",
    signature="",
    definition=function(transition, popSize)
    {
      commute_time <- commute_time_undigraph(transition)
      linearizedFst = commute_time / (16*sum(valuesA(popSize))*ncellA(popSize))
      linearizedFst
    }
  )
    setMethod(
      f="Collisionijk",
      signature="",
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
  
  
  setMethod(
    f="repnDispMutFunction",
    signature="",
    definition=function(geneticData, dimGeneticData, mutationRate, transitionmatrice)
    {
      locusCols = grep("Locus", colnames(geneticData))
      for (individual in 1:dimGeneticData[1]){
        mothercell = sample(nCell, 1,transitionmatrice[geneticData[individual,"Cell_numbers"],])
      geneticline = sample(which(geneticData[,"Cell_numbers"]==mothercell),1)
      geneticData[individual,locusCols] = geneticData[geneticline,locusCols]
    }
    step = 2
    liability = as.data.frame(matrix(liability, ncol = length(locusCols), nrow = dimGeneticData[1]))
    geneticData[,locusCols] = geneticData[,locusCols] + ((liability<mu/2)*step - (liability>(1-mu/2))*step) 
    geneticData
    }

  )
  
  
  setMethod(
    f="conquadraticsq",
    signature="",
    definition=function(X,p)
    {
      xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
      xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]  
      yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]  
      res = (yopt-(4*yopt/(xmax-xmin)^2)*(X-(xmin+xmax)/2)^2)*((X>xmin)&(X<xmax))
      res + res * (1-res)
    }
    
    
  )
  
  
  setMethod(
    f="coalescenceTimeProbDistrib",
    signature="",
    definition=function(Qlist)
    {
      Ndeme <- (2*dim(Qlist[[1]])[1]+1/4)^.5-.5
      cTPD <- array(NA,dim=c(Ndeme,Ndeme,1))
      condition=TRUE
      t=1
      QPower <- Q <- Qlist[["Q"]]
      cTPD[,,t] <- matrix(colSums((matrix(1,nrow(QPower),ncol(QPower))-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme)
      while (condition){
        t=t+1
        QPowerBefore <- QPower
        QPower <- Q%*%QPowerBefore
        cTPD <- abind(cTPD,matrix(colSums((QPowerBefore-QPower)[,c(Qlist[["Qline"]])]),Ndeme,Ndeme),along=3)
        condition=!all(colSums(cTPD,3)>0.9999)
      }
      cTPD
    }
    
    #
    
  )
  
  
  setMethod(
    f="SharedAlleleDistance",
    signature="",
    definition=function(tip_genotype)
    {
      if (all(grep("\\.2",colnames(tip_genotype))==(grep("\\.1",colnames(tip_genotype))+1))) {
        tip_genotype <- TwoCols2OneCol(tip_genotype)
      } 
      Cells <- levels(as.factor(tip_genotype$Cell_numbers))
      shAll_mat <- matrix(NA, nrow=length(Cells),ncol=length(Cells))
      for (cell1 in Cells)
        for(cell2 in Cells)
        {
          alleles1 <- tip_genotype[tip_genotype$Cell_numbers==cell1,]
          alleles2 <- tip_genotype[tip_genotype$Cell_numbers==cell2,]
        }
      g <- tip_genotype[,grep("ocus",colnames(tip_genotype))]
      g1 <- g[,1:ncol(g)/2]
      g2 <- g[,1+1:ncol(g)/2]
    }
    
    #
    #
    
  )
  
  
  setMethod(
    f="SSb",
    signature="",
    definition=function(gDat)
    {
      1+Qwithin_pair(gDat)-2*Qbetween(gDat)
    }
    
  )
  
  
  setMethod(
    f="conquadraticskewedsq",
    signature="",
    definition=function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))))
    {
      Yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]
      Xopt = p[rep("Xopt",dim(X)[1]),colnames(X)]
      Xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]
      Xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
      alpha <- -log(2)/log((Xopt-Xmin)/(Xmax-Xmin))
      Xprime<- ((X-Xmin)/(Xmax-Xmin))^alpha*(Xmax-Xmin)+Xmin
      y <- (Yopt-4*Yopt/((Xmax-Xmin)^2)*(Xprime-(Xmin+Xmax)/2)^2)*((X>=Xmin)&(X<=Xmax))
      y+y*(Yopt-y)
    }
    
    
  )
  
  
  setMethod(
    f="a_matrix",
    signature="",
    definition=function(gDat)
    {
      SSb(gDat)/SSw(gDat)-1/2
    }
    
  )
  
  
  setMethod(
    f="fundamentalMatrixAbsorbingChain",
    signature="",
    definition=function(transientTransitionMatrix)
    { #transtientTransitionMatrix
      Ndeme <- (2*dim(transientTransitionMatrix)[1]+1/4)^.5-.5
      solve(diag(Ndeme*(Ndeme+1)/2)-transientTransitionMatrix)
    }
    
    #
    #
    
  )
  
  
  setMethod(
    f="transitionMatrixForward",
    signature="",
    definition=function(r,K, migration, meth="non_overlap")
    {
      rs = matrix(r,nrow=length(r),ncol=length(r))
      Ku = t(matrix(K,nrow=length(K),ncol=length(K)))
      leave = migration*(1+rs)*t(Ku); leave = leave - diag(leave)
      switch (meth,
              non_overlap = migration * rs * Ku / colSums(rs * t(Ku) * migration),
              overlap = migration * (1+rs) * Ku / (colSums((1+rs) * t(Ku) * migration - t(leave)))
      )
    }  
    
    #
    
  )
  
  
  setMethod(
    f="Qwithin_pop",
    signature="",
    definition=function(gDat)
    {
      gDat <- gDat[,grep("ocus",colnames(gDat))]
      if ((dim(gDat)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat)))) { #gDat is 1_col_diploid
      gDat <- OneCol2TwoCols(gDat)}
    matrix_pair = gDat[,grep(".2", colnames(gDat), fixed = T)]
    matrix_impair = gDat[,grep(".1", colnames(gDat), fixed = T)]
    Qw <- (matrix_pair == matrix_impair)
    Qw = mean(Qw)
    Qw
    }

  )
  
  
  setMethod(
    f="genetDistStepping",
    signature="",
    definition=function(migration,popSize,mutation_rate)
    {
      commute_time <- commute_time_undigraph(migration)
      genetDist = (commute_time / 4 + 2*sum(valuesA(popSize)) )* mutation_rate
      genetDist
    }
    
    
    
  )
  
  
  setMethod(
    f="TwoCols2OneCol_",
    signature="",
    definition=function(gDat)
    {
      notlocusonly <- (length((grep("ocus",colnames(gDat))))<dim(gDat)[2]) 
      if (notlocusonly) {
        coords <- gDat[rep(1:dim(gDat)[2],each=2),-grep("ocus",colnames(gDat))]
        coords <- data.frame(x=rep(gDat[,"x"],each=2),y=rep(gDat[,"y"],each=2),Cell_numbers=rep(gDat[,"Cell_numbers"],each=2))
        rownames(coords)<-paste(paste(rep(rownames(gDat),each=2),each=1:2,sep="."))
      } 
      gDat <- gDat[,grep("ocus",colnames(gDat))]
      pair <- (1:(dim(gDat)[2]/2))*2;impair=pair-1
      first <- gDat[,impair];colnames(first) <- unlist(strsplit(colnames(first),"\\."))[impair]
      second <- gDat[,pair];colnames(second)<-colnames(first)
      Onecol<-rbind(first,second)
      neworder <- rep(c(0,dim(gDat)[1]),dim(gDat)[1])+rep(1:dim(gDat)[1],each=2)
      Onecol <- Onecol[neworder,]
      row.names(Onecol) <- paste(paste(rep(rownames(gDat),each=2),each=1:2,sep="."))
      if (notlocusonly) (cbind(coords,Onecol)) else Onecol
    }
    
  )
  
  
  setMethod(
    f="plot_genealogy",
    signature="",
    definition=function(genealogy,file=NA)
    {
      if (!is.na(file)) pdf(file)
      par(mfrow=c(2,1),oma=c(0,0,0,4),xpd=TRUE)
      plot(coalescent_2_phylog(genealogy),direction="downward")
      if (!is.na(file)) dev.off()
      
    }
    
  )
  
  
  setMethod(
    f="genetDistUndigraphForNLM",
    signature="",
    definition=function(parameters)
    {
      migration <- migrationMatrix(rasterStack = environmentalData,shapeDisp = prior$dispersion$model,pDisp = parameters["dispersion"])
      rasK=environmentalData;rasK[cellNumA(rasK)]= as.matrix(ReactNorm(valuesA(environmentalData),p=parameters["K"],shapes=prior$K$model)[,"Y"]);rasK=rasK[[1]]
      transition <- transitionMatrixBackward(r = parameters["R"],K = valuesA(rasK), migration)
      genetDist_Undigraph <- genetDistUndigraph(transition,popSize=rasK,mutation_rate = parameters["mutation_rate"])
      sum((as.dist(genetDist)-as.dist(genetDist_Undigraph))^2+
          (diag(genetDist)-diag(genetist_Undigraph))^2)/
        (nrow(genetDist)*(nrow(genetDist)+1)/2)
    }
    
  )
  
  
  setMethod(
    f="a_value_matrix",
    signature="",
    definition=function(gDat,rasterStack)
    {
      Cell_numbers = gDat$Cell_numbers
      cell_numbers_levels <- levels(as.factor(gDat$Cell_numbers))
      vecteur_a_value <- a_value_ind(gDat)
      gDat <- gDat[,grepl("\\.",colnames(gDat))]
      if ((dim(gDat)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat)))) {
      cell_numbers_of_ind <- Cell_numbers[(1:(dim(gDat)[1]/2))*2-1]} else {
        cell_numbers_of_ind <- Cell_numbers
      }
    dimnames(vecteur_a_value)<-list(cell_numbers_of_ind,cell_numbers_of_ind)
    agg <- aggregate(vecteur_a_value,by=list(cell_numbers_of_ind),FUN="mean")[,-1]
    agg2 <- t(aggregate(t(agg),by=list(cell_numbers_of_ind),FUN="mean"))[-1,]
    row.names(agg2) <- as.numeric(cell_numbers_levels)
    colnames(agg2) <- as.numeric(cell_numbers_levels)
    agg3 <- matrix(NA, nrow=ncellA(rasterStack),ncol=ncellA(rasterStack))
    agg3[as.numeric(row.names(agg2)),as.numeric(colnames(agg2))] <- agg2
    agg3
    }

  )
  
  
  setMethod(
    f="absorbingTime",
    signature="",
    definition=function(fundamentalMatrixAbsChain,Qline)
    {
      Ndeme <- (2*dim(fundamentalMatrixAbsChain)[1]+1/4)^.5-.5
      absTimesVector <- fundamentalMatrixAbsChain%*%rep(1,nrow(fundamentalMatrixAbsChain))
      absTimesMatrix <- matrix(absTimesVector[Qline],nrow=Ndeme)
      absTimesMatrix
    }
    
    #
    #
    
  )
  
  
  setMethod(
    f="Genetic_Dist",
    signature="",
    definition=function(gDat,method="Goldstein")
    {
      switch(method,
             Goldstein=deltaMu(tip_genotype = gDat,stepvalue = 1)  
      )
    }
    
    #
    #
    
  )
  
  
  setMethod(
    f="enveloppe",
    signature="",
    definition=function(X,p)
    {
      p[rep("Yopt",dim(X)[1]),colnames(X)]*((X>p[rep("Xmin",dim(X)[1]),])&(X<p[rep("Xmax",dim(X)[1]),colnames(X)]))
    }
    
    
  )
  
  
  setMethod(
    f="conquadratic",
    signature="",
    definition=function(X,p)
    {
      xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
      xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]  
      yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]  
      (yopt-(4*yopt/(xmax-xmin)^2)*(X-(xmin+xmax)/2)^2)*((X>xmin)&(X<xmax))
    }
    
    
  )
  
  
  setMethod(
    f="genetDistDigraph",
    signature="",
    definition=function(transition,popSize,mutation_rate,method="Goldstein95")
    {
      H <- hitting_time_digraph(transition)
      dim2 <- dim(H);dim2[[3]]=2
      H2 <- array(c(H,t(H)),dim=dim2)
      MinH <- apply(H2,c(1,2),min)
      genetic_dist = (MinH/2 + 2*sum(valuesA(popSize)))* mutation_rate
      genetic_dist
    }
    
  )
  
  
   
  setMethod(
    f="TwoCols2OneCol",
    signature="",
    definition=function(tip_genotype)
    {
      nonGenotypeData <- tip_genotype[,!grepl("ocus",colnames(tip_genotype))]
      
      matrix_pair = tip_genotype[,grep("\\.2", colnames(tip_genotype))]
      colnames(matrix_pair) <- grep("ocus",unlist(strsplit(colnames(matrix_pair),"\\.")),value=TRUE)
      matrix_pair <- cbind(nonGenotypeData,matrix_pair)
      rownames(matrix_pair)<-paste(rownames(matrix_pair),".1",sep="")
      matrix_impair = tip_genotype[,grep("\\.1", colnames(tip_genotype))]
      colnames(matrix_impair) <- grep("ocus",unlist(strsplit(colnames(matrix_impair),"\\.")),value=TRUE)  
      matrix_impair <- cbind(nonGenotypeData,matrix_impair)
      rownames(matrix_impair)<-paste(rownames(matrix_impair),".2",sep="")  
      tip_genotype <- rbind(matrix_pair,matrix_impair)[rep(1:nrow(matrix_pair),each=2)+rep(c(nrow(matrix_pair),0),nrow(matrix_pair)),]
      tip_genotype
    }
    
  )
  
  
  setMethod(
    f="FstatsRas",
    signature="",
    definition=function(gDat,by="cell",all=TRUE,cells=NULL)
    {
      Absent <- !cells%in%gDat$Cell_numbers
      if ((all==TRUE)&(any(Absent))){
        AbsentCells <- (cells)[Absent]
        gDat[dim(gDat)[1]+1:length(AbsentCells),"Cell_numbers"] <-AbsentCells
      }
      for (i in levels(as.factor(gDat$Cell_numbers))){
        gDat[gDat$Cell_numbers==as.numeric(i),"newCelN"] <- which(levels (as.factor(gDat$Cell_numbers))==i)
      }
      gDat2 <- gDat[,grep("ocus",colnames(gDat))]
      if ((dim(gDat2)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat2)))) {
      gDat <- OneCol2TwoCols(gDat)}
    genotypes <- gDat[,grep("ocus",colnames(gDat))]
    cells <- as.integer(levels(as.factor(gDat$Cell_numbers)))
    Fstat(genotypes,npop=length(cells),pop.mbrship=gDat$newCelN,ploidy=2)
    }

#
#

  )


setMethod(
  f="coalescent_2_newick",
  signature="",
  definition=function(coalescent)
  {
    tree=paste(" ",coalescent[[length(coalescent)]]$new_node," ",sep="")
    for (i in length(coalescent):1)
    {
      Time = coalescent[[i]]$time
      coalesc <- as.character(coalescent[[i]]$coalescing)
      tree <- str_replace(tree,paste(" ",as.character(coalescent[[i]]$new_node)," ",sep=""),paste(" ( ",paste(" ",coalesc," :",coalescent[[i]]$br_length,collapse=" ,",sep=""),") ",sep=""))
    }
    tree <- gsub(" ","",paste(tree,";",sep=""))
    tree
  }
  
  #
  
)


setMethod(
  f="expect_linearizedFst_undigraph",
  signature="",
  definition=function(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA)
  {
    if (is.na(Cell_numbers[1])) Cell_numbers = cellNumA(rasterStack)
    K = ReactNorm(valuesA(rasterStack), pK, shapesK)[,"Y"]
    r = ReactNorm(valuesA(rasterStack), pR, shapesR)[,"Y"]
    rasK=rasterStack;rasK[cellNumA(rasK)]= as.matrix(ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"])
    matrice_migration = migrationMatrix(rasK,shapeDisp, pDisp)
    matrice_transition = transitionMatrixBackward(r,K, matrice_migration)
    linearizedFst = linearizedFstUndigraph(matrice_transition, popSize=rasK)[Cell_numbers,Cell_numbers]
    linearizedFst
  }
  
)


setMethod(
  f="coalist_2_coaltable",
  signature="",
  definition=function(coalist,locusnames,initial_locus_value)
  {
    coaldf <- data.frame(Reduce(rbind,coalist))
    coaltable <- coaldf[rep(1:dim(coaldf)[1],unlist(lapply(coaldf$coalescing, length))),]
    coaltable[,c("coalescing","br_length",locusnames)] <- unlist(coaldf[,c("coalescing","br_length",locusnames)])
    coaltable$genotype <- NA
    coaltable[dim(coaltable)[1]+1,"genotype"] <- initial_locus_value
    coaltable[dim(coaltable)[1],"coalescing"] <- dim(coaltable)[1]
    coaltable$mut <- coaltable[,grep("Locus",colnames(coaltable),value=TRUE)]
    coaltable <- coaltable[,-which(colnames(coaltable)%in%grep("Locus",colnames(coaltable),value=TRUE))]
  }
  
  #
  
)


setMethod(
  f="nbpDisp",
  signature="",
  definition=function(shapeDisp)
  {
    (switch(shapeDisp,
            fat_tail1 = 2,
            gaussian = 1,
            exponential = 1,
            contiguous = 1,
            island = 1,
            fat_tail2 = 2))
  }
  
)


setMethod(
  f="ReactNorm",
  signature="",
  definition=function(X,p,shapes)
  {
    if (class(p)=="numeric") {p=t(as.matrix(p));rownames(p)="a"}
    if (class(X)=="numeric") {X=data.frame(X=X);colnames(X)=colnames(p)}
    Y=X
    if (!all(colnames(p)%in%names(shapes))) {stop ("variable names do not correspond between parameters 'p' ans 'shape'")}
    for (shape in as.character(levels(as.factor(shapes))))
    {
      variables = colnames(p)[which(shapes==shape)]
      Y[,variables]=switch(shape,
                           constant=p,
                           proportional = proportional(subset(X,select=variables),p),
                           linear = linear(subset(X,select=variables),p),
                           enveloppe=enveloppe(subset(X,select=variables),p),
                           envelin=envelinear(subset(X,select=variables),p),
                           envloglin=envelinear(subset(X,select=variables),p,log=TRUE),
                           loG = log(subset(X,select=variables)),
                           conquadratic=conquadratic(subset(X,select=variables),p),
                           conquadraticskewed=conquadraticskewed(subset(X,select=variables),p),
                           conquadraticsq=conquadraticsq(subset(X,select=variables),p),
                           conquadraticskewedsq=conquadraticskewedsq(subset(X,select=variables),p)
      )
    }
    Y
  }
  
  
)


setMethod(
  f="migrationMatrix",
  signature="",
  definition=function(rasterStack,shapeDisp, pDisp)
  {
    distMat<-distanceMatrix(rasterStack)
    Ndim = 1+all(ncell(rasterStack)!=dim(rasterStack)[1:2])
    migration = apply(distMat, c(1,2), 
                      function(x)(switch(shapeDisp,
                                         gaussian = (dnorm(x, mean = 0, sd = pDisp[1], log = FALSE)),
                                         exponential = (dexp(x, rate = 1/pDisp[1], log = FALSE)),
                                         contiguous = (x==0)*(1-pDisp[1])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp[1]/(2*Ndim)),
                                         contiguous8 = (x==0)*(1-pDisp[1])+((x>0)-(x>2*res(rasterStack)[1]))*(pDisp[1]/(4*Ndim)),
                                         island = (x==0)*(1-pDisp[1])+(x>0)*(pDisp[1]),
                                         fat_tail2 = x^pDisp[2]*exp(-2*x/(pDisp[1]^0.5)),
                                         contiguous_long_dist_mixt = pDisp["plongdist"]/ncellA(rasterStack)+(x==0)*(1-pDisp["pcontiguous"]-pDisp["plongdist"])+((x>0)-(x>1.4*res(rasterStack)[1]))*(pDisp["pcontiguous"]/2),
                                         gaussian_long_dist_mixt = pDisp[2]/ncellA(rasterStack) + (dnorm(x, mean = 0, sd = pDisp[1], log = FALSE))
                      )))
    return(migration)
  }
  
)


setMethod(
  f="absorbingTransitionComplete",
  signature="",
  definition=function(transition,N)
  {
    N[N<1]=1
    Ndeme <- dim(transition)[1]
    Nhetero <- Ndeme*(Ndeme-1)/2
    Qline <- matrix(NA,Ndeme,Ndeme) # this matrix of the same dimension as the transition matrix
    HoHe <- matrix(NA,Ndeme,Nhetero)
    HeHo <- matrix(NA,Nhetero,Ndeme)
    CoHe <- matrix(0,Ndeme,Nhetero)
    HeCo <- matrix(NA,Ndeme,Nhetero)
    HoCo <- matrix(transition/(2*matrix(N,nrow=Ndeme,ncol=Ndeme, byrow=TRUE)),nrow=Ndeme,ncol=Ndeme)
    CoHo <- matrix(0,nrow=Ndeme,ncol=Ndeme)
    HeHe <- matrix(NA,Nhetero,Nhetero)
    CoCo <- transition
    HoHo <- transition*transition
    ij=0;kl=0
    Check=TRUE
    for (i in 2:Ndeme){
      for(j in 1:(i-1)) {
        ij=ij+1
        Qline[i,j] <- Qline[j,i] <- ij # this matrix aims to find where the pair of demes {i,j}
        kl=0
        for (k in 2:Ndeme){
          for (l in 1:(k-1)){
            kl=kl+1
            HeHe[ij,kl] <- transition[i,k]*transition[j,l]
          }
        }
        for (k in 1:Ndeme){
          HeHo[ij,k] <- transition[i,k]*transition[j,k]*(1-1/(2*N[k]))
          HeCo[ij,k] <- transition[i,k]*transition[j,k]/(2*N[k])
          HoHe[k,ij] <- transition[k,i]*transition[k,j]
        }
      }
    }
    
    kl=0
    Check=TRUE
    for (i in 1:Ndeme){ # only homodemic sources are considered (i=j)
      Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
      for (k in 2:Ndeme){
        for (l in 1:(k-1)){ # only heterodemic targets
          kl=kl+1
    }
  }
  kl=0
  }
QheteroHomo <- matrix(NA,Nhetero,Ndeme)
HeCo <- matrix(NA,nrow=Nhetero,ncol=Ndeme)
ij=0
for (i in 2:Ndeme){
  for (j in 1:(i-1)){ # only heterodemic sources
    ij=ij+1
    for(k in 1:Ndeme){ # only homodemic targets
      QheteroHomo[ij,k] <- transition[i,k]*transition[j,k]*(1-1/(2*N[k]))
      HeCo[ij,k] <- transition[i,k]*transition[j,k]/(2*N[k])
    }
  }
}
QhomoHomo <- transition*transition*matrix(1-1/(2*N),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
Q <- cbind(rbind(QheteroHetero,QhomoHetero),rbind(QheteroHomo,QhomoHomo))
CoHeCoHo <- matrix(0,nrow=Ndeme,ncol=Nhetero+Ndeme)
ij=0
list(Q=Q,Qline=Qline)
}

#
#

)


setMethod(
  f="samplePrior",
  signature="",
  definition=function(prior,method="random")
  {
    p <- NA
    if (method=="random") {
      for (i in names(prior)){
        p[i] = switch(prior[[i]]$distribution,
                      uniform = runif(1,min=prior[[i]]$p["min"],max=prior[[i]]$p["max"]),
                      loguniform = exp(runif(1,min=log(prior[[i]]$p["min"]),max=log(prior[[i]]$p["max"]))),
                      fixed = prior[[i]]$p)
      }
    }
    if (method=="mean") {
      for (i in names(prior)){
        p[i] = switch(prior[[i]]$distribution,
                      uniform = mean(c(prior[[i]]$p["min"],prior[[i]]$p["max"])),
                      loguniform = exp(mean(c(min=log(prior[[i]]$p["min"]),max=log(prior[[i]]$p["max"])))),
                      fixed = prior[[i]]$p)
      }
    }
    p[-1]
  }
  
)


setMethod(
  f="ocur",
  signature="",
  definition=function(cells=c(1,2),transitionmatrice)
  {
    tmpcell1 <- cells[1];tmpcell2 <- cells[2];t=1
    while (tmpcell1 != tmpcell2)
    {
      tmpcell1 =  sample(dim(transitionmatrice)[2],size=1,prob=c(transitionmatrice[tmpcell1,]))
      tmpcell2 =  sample(dim(transitionmatrice)[2],size=1,prob=c(transitionmatrice[tmpcell2,]))
      t=t+1
    }
    t
  }
  
  #
  #
)


setMethod(
  f="degree2km",
  signature="",
  definition=function(rasterStack)
  {
    
    dist_degree <- acos(sin(x_origin)*sin(x_destination)+cos(x_origin)*cos(x_destination)*cos(y_origin-y_destination))
    dist_km = dist_degree * 111.32
    dist_km
  }
  
)


setMethod(
  f="checkTimeInterval",
  signature="",
  definition=function(transition = result$transitionmatrice,rasK = rasK,threshold = 1E-5)
  {
    meanCoalTimes <- list()
    for (time_interval in c("1","4","10","15","20","40","100")){
      meanCoalTimes[[time_interval]] <- coalescence_prob_time_distribution_matrix(transition = transition,max_time_interval = as.numeric(time_interval),rasK = rasK,threshold = 1E-5)$exp_times
    }
    meanCoalTimes
  }
  
  #
  #
  
)


setMethod(
  f="conquadraticskewed",
  signature="",
  definition=function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))))
  {
    Yopt = p[rep("Yopt",dim(X)[1]),colnames(X)]
    Xopt = p[rep("Xopt",dim(X)[1]),colnames(X)]
    Xmin = p[rep("Xmin",dim(X)[1]),colnames(X)]
    Xmax = p[rep("Xmax",dim(X)[1]),colnames(X)]
    alpha <- -log(2)/log((Xopt-Xmin)/(Xmax-Xmin))
    Xprime<- ((X-Xmin)/(Xmax-Xmin))^alpha*(Xmax-Xmin)+Xmin
    y <- (Yopt-4*Yopt/((Xmax-Xmin)^2)*(Xprime-(Xmin+Xmax)/2)^2)*((X>=Xmin)&(X<=Xmax))
    y[X<Xmin] <-0
    y
  }
  
)


setMethod(
  f="expect_linearizedFst_digraph",
  signature="",
  definition=function(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA)
  {
    if (is.na(Cell_numbers[1])) Cell_numbers = cellNumA(rasterStack)
    K = ReactNorm(valuesA(rasterStack), pK, shapesK)[,"Y"]
    r = ReactNorm(valuesA(rasterStack), pR, shapesR)[,"Y"]
    rasK=rasterStack;rasK[cellNumA]= as.matrix(ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"])
    matrice_migration = migrationMatrix(rasK,shapeDisp, pDisp)
    matrice_transition = transitionMatrixBackward(r,K, matrice_migration)
    linearizedFst_att = linearizedFstDigraph(matrice_transition,rasK)[Cell_numbers,Cell_numbers]
    linearizedFst_att
  }
  
)


setMethod(
  f="deltaMu",
  signature="",
  definition=function(tip_genotype,stepvalue=2)
  {
    g <- aggregate(tip_genotype,by=list(tip_genotype$Cell_numbers),FUN="mean",na.action="na.omit")
    Cells <-g$Cell_numbers  
    g <- g[,grep("ocus",colnames(tip_genotype))]
    gd <- as.matrix(dist(g/stepvalue)^2/dim(g)[2])
    dimnames(gd) <- list(Cells,Cells)
    gd
  }
  
  #
  #
  
)


setMethod(
  f="Qwithin_pair_2",
  signature="",
  definition=function(gDat)
  {
    gDatG <- gDat[,grep("ocus",colnames(gDat))]
    gDat3D <- OneCol23Dims(gDatG)
    (t(array(diag(id_mat_prob(gDat3D[,,1],gDat3D[,,2])),dim=c(dim(gDat3D)[1],dim(gDat3D)[1])))+array(diag(id_mat_prob(gDat3D[,,1],gDat3D[,,2])),dim=c(dim(gDat3D)[1],dim(gDat3D)[1])))/2
  }
  
)


setMethod(
  f="envelin0",
  signature="",
  definition=function(X,p,log=FALSE)
  {
    Yxmin = p[rep("Yxmax",dim(X)[1]),colnames(X)]
    Yxmax = p[rep("Yxmin",dim(X)[1]),colnames(X)]
    Xmin=  p[rep("Xmin",dim(X)[1]),colnames(X)]
    Xmax=  p[rep("Xmax",dim(X)[1]),colnames(X)]
    a = (Yxmin - Yxmax) / (Xmin - Xmax)
    b = Yxmin - Xmin * a
    if (log) {log(a*X+b)*((X>Xmin)&(X<Xmax))} else {(a*X+b)*((X>Xmin)&(X<Xmax))}
  }
  
)


setMethod(
  f="transitionMatrixBackward",
  signature="",
  definition=function(r,K, migration)
  {
    if ((length(r)==1)&(length(K)==1)){transition = r * K * t(migration)}
    if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
    if ((length(r)==1)&(length(K)>1)){transition = r * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
    if ((length(r)>1)&(length(K)==1)){transition = t(matrix(r,nrow=length(r),ncol=length(r))) * K * t(migration)}
    if ((length(r)>1)&(length(K)>1)) {transition = t(matrix(r,nrow=length(r),ncol=length(r))) * t(matrix(K,nrow=length(K),ncol=length(K))) * t(migration)}
    transition
  }
  
)


setMethod(
  f="validation",
  signature="",
  definition=function(donneesEnvironmentObs,pK = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))), pR = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))),shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"),shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"), shapeDisp="fat_tail1", pDisp = c(mean=0.32,shape=1.6), file=NULL,mutationRate,nbLocus, initial_locus_value,nb_generations=5000,indpercell=30)
  {
    nblayers =dim(donneesEnvironmentObs)[3]
    nCell = ncellA(donneesEnvironmentObs)
    Cell_numbers <- cellNumA(donneesEnvironmentObs)
    K <- ReactNorm(valuesA(donneesEnvironmentObs), pK , shapesK)$Y # recuperation de la carrying capacity
    r <- ReactNorm(valuesA(donneesEnvironmentObs), pR , shapesR)$Y # recuperation du growth rate
    migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
    transitionmatrice = transitionMatrixBackward(Npop = K, migration= migrationM)
    geneticObs = simulationGenet(donneesEnvironmentObs,alpha, beta, mutationRate,nbLocus, initial_locus_value,shapeDisp,pDisp,nb_generations=5000,indpercell)
    finalGenetData = geneticObs[[1]]
    land_size <- raster(matrix(K[1,],nrow=dim(donneesEnvironmentObs)[1],ncol=dim(donneesEnvironmentObs)[2]))
    a_value_simul = a_value_matrix(geneticObs[[1]],land_size)
    
    a_value_theory_stepping_stone_model = S/(32*dnorm(res(rasK)[1],0,pDisp[1])*sum(valuesA(rasK))/ncellA(rasK))
    distanceMatrix(donneesEnvironmentObs)/(4*0.05)
    a_value_theory_island_model <- matrix(pDisp[1],nrow=nCell,ncol=nCell)-pDisp[1]*diag(nCell)
    a_value_theory_graph_model <- expect_linearizedFst_undigraph(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers)
    if (!is.null(file)){jpeg(file)}
    par(mfrow=c(2,3))
    plot(land_size,main="Population sizes")
    plot(raster(transitionmatrice),main="Transition matrix")
    plot(raster(a_value_simul),main="Simul genet differ")
    plot(raster(a_value_theory_stepping_stone_model),main="Expected stepping stone")
    plot(raster(a_value_theory_island_model),main="Expected island")
    plot(raster(a_value_theory_graph_model),main="Expected graph model")
    if (!genetDistis.null(file)){dev.off()}
    list(land_size,transitionmatrice,a_value_simul,a_value_theory_stepping_stone_model,a_value_theory_island_model,a_value_theory_graph_model,finalGenetData)
  }
)
  
setMethod(
  f="validation_with_coalescent",
  signature="",
  definition= function(rasterStack=rasterStack,pK = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"), c("BIO1","BIO12"))), pR = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1), nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))),shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"),shapesR=c(BIO1="constant",BIO12="constant"), shapeDisp="gaussian", pDisp = .5, filen=NULL,nbLocus=60, initial_locus_value=200,stepvalue=2,mutation_model="stepwise",mutation_rate=1E-2,nind="Ne", stat="a_value",Option = "sample_1col_diploid")
  {
    rasK=rasterStack;rasK[cellNumA(rasK)]= as.matrix(ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"]);rasK=rasK[[1]]
    geneticData = CreateGenetArray(rasK, nbLocus,initial_locus_value,Option=Option,nind=nind)
    coal_genet <- simul_coal_genet(geneticData,rasterStack,pK=pK,pR=pR,shapesK=shapesK,shapesR=shapesR,
                                   shapeDisp=shapeDisp,pDisp=pDisp,
                                   mutation_rate=mutation_rate,
                                   initial_locus_value=initial_locus_value,
                                   mutation_model=mutation_model,
                                   stepvalue=stepvalue,locusnames=NA)
    a_valueSimul <- a_value_matrix(coal_genet$tip_genotype,rasK) 
    FstSimul <- FstatsRas(coal_genet$tip_genotype,by="cell",all=TRUE,cells=cellNumA(rasterStack))$Fst
    genetDistSimul <- genetDist(tip_genotype = coal_genet$tip_genotype,method = "deltaMu",stepvalue = 2,byCell = TRUE)
    if (shapeDisp %in% c("contiguous","gaussian")){
      stepping_migration_rate <- switch(shapeDisp,
                                        contiguous=pDisp,
                                        gaussian=2*dnorm(res(rasK)[1],0,pDisp))
      genetDistStepGM <- genetDistStepping(migrationMatrix(rasK,shapeDisp,pDisp),popSize = rasK,mutation_rate = mutation_rate)
      genetDistIslandGM <- matrix(1/(1+4*stepping_migration_rate*mean(valuesA(rasK)))/(1-1/(1+4*stepping_migration_rate*mean(valuesA(rasK)))),nrow=ncellA(rasK),ncol=ncellA(rasK))
      diag(genetDistIslandGM) <- 1/(1+4*(1-stepping_migration_rate*mean(valuesA(rasK))))
      if (any(dim(rasterStack) == ncell(rasK))) {
      linearizedFst_theory_stepping_stone_model = distanceMatrix(rasK)/(res(rasK)[1])
      genetDistStepM <- 1/(stepping_migration_rate/2*distanceMatrix(rasK)/res(rasK))^2
      linearizedFst_theory_stepping_stone_model = log(distanceMatrix(rasterStack))/(res(rasterStack)[1])
      genetDistStepM <- log(1/(stepping_migration_rate/4*distanceMatrix(rasK)/res(rasK))^2)
    }
  } else {
    linearizedFst_theory_stepping_stone_model=matrix(NA,nrow=dim(transition)[1],ncol=dim(transition)[2])
    genetDistStepM=NA 
    genetDistStepMG=NA
  }
  linearizedFst_theory_island_model <- matrix(pDisp[1],nrow=ncellA(rasK),ncol=ncellA(rasK))-pDisp[1]*diag(ncellA(rasK))
  linearizedFst_theory_undigraph_model <- expect_linearizedFst_undigraph(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA)
  genetDist_Undigraph <- genetDistUndigraph(coal_genet$transition$backw,popSize=rasK,mutation_rate = mutation_rate)
  genetDist_Digraph <- genetDistDigraph(coal_genet$transition$backw,popSize=rasK,mutation_rate = mutation_rate)
  linearizedFst_theory_digraph_model <- expect_linearizedFst_digraph(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA)
  coalescence_times <- coalescence_prob_time_distribution_matrix(transition = coal_genet$transitionAndK$backw,max_time_interval = 4,rasK = rasK,threshold = 1E-5)
  coalescence_distance <- coalescence_times[["exp_times"]]*mutation_rate
  hitting_time <- hitting_time_digraph(coal_genet$transitionAndK$backw)
  sampleSize<- rasK
  cellsSampled <- as.numeric(levels(as.factor(geneticData$Cell_numbers)))
  for (i in cellNumA(rasK)) {values(sampleSize)[i] <- sum(coal_genet$tip_genotype$Cell_numbers==i)}
  genetDist_absorbingMethod <- genetDistAbsorbingMethod(coal_genet$transitionAndK$backw,valuesA(rasK),mutation_rate) 
  rasK<-as.matrix(rasK)
  
  result = list(coalescents_simulated=coal_genet$coalescent,
                tip_genotypes=coal_genet$tip_genotype,
                rasK=rasK,
                popSize = as.matrix(rasK),
                sampleSize=as.matrix(sampleSize),
                cellsSampled=cellsSampled,
                transitionmatrice=coal_genet$transitionAndK$backw,
                a_valueSimul=a_valueSimul,
                Fst_Simul=FstSimul,
                linearizedFstSimul=FstSimul/(1-FstSimul),
                geneticDistSimul = genetDistSimul,
                linearizedFst_theory_undigraph_model=linearizedFst_theory_undigraph_model,
                linearizedFst_theory_digraph_model=linearizedFst_theory_digraph_model,
                genetDistUndigraph=  genetDist_Undigraph,
                genetDistDigraph=  genetDist_Digraph,
                hitting_time=hitting_time,
                coalescence_times_distribution=coalescence_times[[1]],
                coalescence_times_mean=coalescence_times[[2]],
                coalModelMutDist=coalescence_distance,
                genetDist_absorbingMethod=genetDist_absorbingMethod,
                parameters = list(rasterStack=rasterStack,pK=pK,pR =pR,
                                  shapesK=shapesK,shapesR=shapesR, shapeDisp=shapeDisp, 
                                  pDisp = pDisp, filen=filen,mutation_rate=mutation_rate, nbLocus=nbLocus, 
                                  initial_locus_value=initial_locus_value,stepvalue=stepvalue,
                                  mutation_model=mutation_model,mutation_rate=mutation_rate,
                                  nind=nind,stat=stat,Option=Option)
  )
  if (shapeDisp %in% c("contiguous","gaussian")){result[["linearizedFst_theory_stepping_stone_model"]]=linearizedFst_theory_stepping_stone_model
  result[["genetDistStepping"]]=genetDistStepM
  result[["genetDistGraphStepping"]]=genetDistStepGM
  result[["genetDistIsland"]]=genetDistIslandGM
  }
  result
  }

#
#
#
#

)


setMethod(
  f="transitionMatrice",
  signature="",
  definition=function(rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp)
  {
    K = ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"]
    r = ReactNorm(valuesA(rasterStack),pR,shapesR)[,"Y"]
    migrationM <- migrationMatrix(rasterStack,shapeDisp, pDisp)
    transitionmatrice = transitionMatrixBackward(r, K, migration= migrationM);
    transition_forward = transitionMatrixForward(r, K, migration= migrationM)
    list(backw=transitionmatrice,forw=transition_forward,K=K)  
  }
  
  
)


setMethod(
  f="genetDistUndigraph",
  signature="",
  definition=function(transition,popSize,mutation_rate)
  {
    commute_time <- commute_time_undigraph(transition)
    genetDist =  (commute_time/4+2*sum(valuesA(popSize))) * mutation_rate
    genetDist
  }
  
)


setMethod(
  f="Qwithin_pair",
  signature="",
  definition=function(gDat)
  {
    gDat <- gDat[,grep("ocus",colnames(gDat))]
    if ((dim(gDat)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat)))) { #gDat is 1_col_diploid
    gDat <- OneCol2TwoCols(gDat)}
  matrix_pair = gDat[,grep(".2", colnames(gDat), fixed = T)]
  matrix_impair = gDat[,grep(".1", colnames(gDat), fixed = T)]
  Qw <- (matrix_pair == matrix_impair)
  Qw = rowMeans(Qw) # vector of probability of Qw for each individual
  Qw = matrix(Qw, ncol = dim(gDat)[1], nrow = dim(gDat)[1])
  Qw = (Qw+t(Qw))/2
  Qw
  }

)


setMethod(
  f="a_value_ind",
  signature="",
  definition=function(gDat)
  {
    matrixQb = (1-Qwithin_pair(gDat)+2*(Qwithin_pair(gDat)-Qbetween(gDat)))
    matrixQw = 2*(1-Qwithin_pop(gDat))
    vecteur_a_value = matrixQb/matrixQw-1/2
    vecteur_a_value[is.na(vecteur_a_value)] <-0
    vecteur_a_value
  }
  
  
  #
  #
  
  niche_model <- function(rasterStack,
                          abundanceData,
                          priorsK = list(min=matrix(c(10,200,11,2,0,0,250,1500,1500,2,0,0),
                                                    nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),
                                                                                c("BIO1","BIO12")))
                                         , max=matrix(c(200,400,300,30,0,0,300,2000,1500,30,0,0),
                                                      nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),
                                                                                  c("BIO1","BIO12")))),
                          shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"))
  {
    rasK=rasterStack;valuesA(rasK)= as.matrix(ReactNorm(valuesA(rasterStack),pK,shapesK)[,"Y"]);rasK=rasK[[1]]
    abundanceData <- cbind(abundancedata,extract(rasK,spatialPoints(abundanceData[,c("x","y")])))
  }
  
  
  
)


setMethod(
  f="Qbetween",
  signature="",
  definition=function(gDat)
  {
    colNames  <- colnames(gDat)
    gDat <- gDat[,grep("ocus",colnames(gDat))]
    if ((dim(gDat)[1] %% 2 == 0) & !all(grepl("\\.",colnames(gDat)))) {
    gDat <- OneCol2TwoCols(gDat)}
    Qb = matrix(ncol = dim(gDat)[1], nrow = dim(gDat)[1]) #initialization of Qb as a matrix
  A=as.matrix(gDat[,grep("ocus",colnames(gDat),fixed=T)])
  A3 = aperm(array(A,dim=c(dim(A)[1],dim(A)[2],dim(A)[1])), c(1,3,2)) # permutation des colonnes et des etages
  B3 = aperm(A3, c(2,1,3)) # transposee de A3
  moy1 = colMeans(aperm(A3 == B3), dims = 1, na.rm = T)
  
  l= 1:dim(A)[2]
  Aprime= A[,c(matrix(c(l[2*floor(1:(length(l)/2))],l[2*floor(1:(length(l)/2))-1]), nrow= 2, ncol = length(l)/2, byrow = T))] # Permutation des colonnes
  Aprime3 = aperm(array(Aprime,dim=c(dim(A)[1],dim(A)[2],dim(A)[1])), c(1,3,2))
  moy2 = colMeans(aperm(Aprime3 == B3), dims = 1, na.rm = T) # calcul moy dist pour les individus avec loci permutés Aprime3 et la transposée B3
  Qb =(moy1 + moy2)/2
  Qb
  }

)


setMethod(
  f="add_br_length_and_mutation",
  signature="",
  definition=function(coalescent,mutation_rate,initial_locus_value,allelenames)
  {
    tips = NULL
    internals = NULL
    nodes = NULL
    times = NULL
    {
      nodes = append(nodes,c(coalescent[[i]]$coalescing,coalescent[[i]]$new_node))
      internals = append(internals,coalescent[[i]]$new_node)
      times = append(times,coalescent[[i]]$time)
    }
    nodes = as.numeric(levels(as.factor(c(nodes,internals))));nodes = nodes[order(nodes)]
    tips = nodes[!((nodes)%in%(internals))]
    {
      {
        if (coalescing %in% tips) {coalescent[[i]]$br_length <- append(coalescent[[i]]$br_length,coalescent[[i]]$time)
        } else {
          coalescent[[i]]$br_length <- append(coalescent[[i]]$br_length,coalescent[[i]]$time-times[which(internals==coalescing)]) 
        } 
        for (allele in allelenames)
        {
          coalescent[[i]][[allele]] <- rpois(rep(1,length(coalescent[[i]]$br_length)),coalescent[[i]]$br_length*mutation_rate)
        }
      }
    }
    coalescent
  }
  
  
  
  
)


setMethod(
  f="distanceMatrix",
  signature="",
  definition=function(rasterStack)
  {
    coords = xyFromCellA(rasterStack)
    return(distance)
  }
  
)


setMethod(
  f="commute_time_undigraph",
  signature="",
  definition=function(matrice_transition)
  {
    laplacian = laplaceMatrix(matrice_transition)
    mii = matrix(diag, nrow =dim(inverseMP), ncol = dim(inverseMP))
    mjj = t(mii)
    mij = inverseMP
    mji = t(mij)
    commute_time = mii + mjj - mij - mji
    commute_time
  }
  
  
)


setMethod(
  f="formatGeneticData",
  signature="",
  definition=function(gDat,rasK)
  {
    loci<- colnames(gDat)[grep("ocus",colnames(gDat))]
    ProbCellPairs<- array(0,dim=c(dim(gDat)[1]/2,length(loci),2,ncellA(rasK)),dimnames=list(paste("pair",1:(dim(gDat)[1]/2),sep=""),loci,c("gene1","gene2"),cellNumA(rasK)))
    GenetDistPairs <- array(NA,dim=c(dim(gDat)[1]/2,length(loci)),dimnames=list(paste("pair",1:(dim(gDat)[1]/2),sep=""),loci))
    for (locus in loci){
      gene1 <- sample(dim(gDat)[1],dim(Pairs)[1])
      gene2 <- sample((1:dim(gDat)[1])[-gene1])
      ProbCellPairs[,locus,1,gDat[gene1,"Cell_numbers"]] <- 1
      ProbCellPairs[,locus,2,gDat[gene2,"Cell_numbers"]] <- 1
      GenetDistPairs[,locus] <- abs(gDat[gene1,locus]-gDat[gene2,locus])
    }
    list(ProbCellPairs=ProbPairs,GenetDistPairs=GenetDistPairs)
  }
  
)


setMethod(
  f="absorbingTransition",
  signature="",
  definition=function(transition,N)
  {
    N[N<1]=1
    Ndeme <- dim(transition)[1]
    Nhetero <- Ndeme*(Ndeme-1)/2
    Qline <- matrix(NA,Ndeme,Ndeme)
    QheteroHetero <- matrix(NA,Nhetero,Nhetero)
    ij=0;kl=0
    Check=TRUE
    for (i in 1:Ndeme){ #only homodemic sources are considered (i=j)
      Qline[i,i] <- Nhetero+i # the homodemic states are after the  heterodemic states
      # in the lines of Q matrix
      for (k in 2:Ndeme){
        for (l in 1:(k-1)){ # only heterodemic targets
          kl=kl+1
          QhomoHetero[i,kl] <- transition[i,k]*transition[i,l]    # i=j (homodemic sources)
          #          Check <- Check*(P_i_k_j_l[P_ij_kl[ij,kl,"i"],P_ij_kl[ij,kl,"j"],P_ij_kl[ij,kl,"k"],P_ij_kl[ij,kl,"l"]]==Q[ij,kl])
        }
      }
      kl=0
    }
    QheteroHomo <- matrix(NA,Nhetero,Ndeme)
    ij=0
    for (i in 2:Ndeme){
      for (j in 1:(i-1)){ # only heterodemic sources
        ij=ij+1
        for(k in 1:Ndeme){ # only homodemic targets
          QheteroHomo[ij,k] <- transition[i,k]*transition[j,k]*(1-1/(2*N[k]))
          # homodemic targets that have not coalesced
        }
      }
    }
    QhomoHomo <- transition*transition*matrix(1-1/(2*N),nrow=Ndeme,ncol=Ndeme,byrow=TRUE)
    Q <- cbind(rbind(QheteroHetero,QhomoHetero),rbind(QheteroHomo,QhomoHomo))
    list(Q=Q,Qline=Qline)
  }
)


setMethod(
  f="landGenetAnalysis",
  signature="",
  definition=function(genetData,environmentalData,priors)
  {
    p=samplePrior(prior,method="random")
    gDist <- genetDist(genetData,method="deltaMu",stepvalue)
    genetDistUndigraphForNLM(p)
    parameters
  }
  
  
  #
  
)


setMethod(
  f="plot_coal_time_depending_on_time_interval",
  signature="",
  definition=function(meanCoalTimes,filen=NA,corr=TRUE,land=TRUE)
  {
    if (land){
      if (!is.na(filen)) {
        if (grepl(".jpg",filen)) jpeg(paste("land_",filen,sep="")) else if (grepl(".pdf",filen)) pdf(paste("land_",filen,sep=""))
      }
      par(mfrow=c(3,3))
      for (i in 1:length(meanCoalTimes)){
        plot(raster(meanCoalTimes[[i]]),main=paste("jump",names(meanCoalTimes)[i]))
      }
      if (!is.na(filen)) dev.off()
    }
    if (corr){
      df<-as.data.frame(matrix(unlist(meanCoalTimes),nrow=prod(dim(meanCoalTimes[[1]])),ncol=length(meanCoalTimes),dimnames=list(1:prod(dim(meanCoalTimes[[1]])),names(meanCoalTimes))))
      if (!is.na(filen)) {
        if (grepl(".jpg",filen)|grepl(".jpeg",filen)) jpeg(paste("corr_",filen,sep="")) 
        if (grepl(".pdf",filen)) pdf(paste("corr_",filen,sep=""))
      }
      plot(df,cex=.1)    
      if (!is.na(filen)) dev.off()
      if (!is.na(filen)) {
        if (grepl(".jpg",filen)) jpeg(paste("corr_",filen,sep="")) 
        if (grepl(".pdf",filen)) pdf(paste("corr_",filen,sep=""))
      }
      if (!is.na(filen)) write.table(cor(df),file=paste("corr_",strsplit(filen,"\\.")[[1]][1],".csv",sep=""))    
      if (!is.na(filen)) dev.off()
    }
  }
  
  
  #
  #
  
  
  plot_validation <- function(validationdata,what=c("pop_size","sample_size","transition","geneticDistSimul",
                                                    "hitting_time","genetDistUndigraph",
                                                    "genetiDistDigraph","graphCoalTimes","graphCoalModelMutDist",
                                                    "Fsimul","FlinearSimul",
                                                    "linFstepping","genetDistStepping",
                                                    "genetDistIsland","FlinearDigraph",
                                                    "FlinearUndigraph","genetDistGraphStepping",
                                                    "genetDistAbsorbingMC"),
                              whichCell="all",
                              filen=NA)
  {
    refWhat <- c("pop_size","sample_size","transition","geneticDistSimul",
                 "hitting_time","genetDistUndigraph",
                 "genetiDistDigraph","graphCoalTimes","graphCoalModelMutDist",
                 "Fsimul","FlinearSimul",
                 "linFstepping","genetDistStepping","genetDistIsland","FlinearDigraph",
                 "FlinearUndigraph","genetDistGraphStepping","genetDistAbsorbingMC")
    if (is.numeric(what)) what <- refWhat[what]
    par(mar=c(3,3,3,3))
    cells<-as.numeric(levels(as.factor(validationdata$tip_genotypes$Cell_numbers)))
    n1=ceiling(sqrt(length(what)))
    n2=ceiling(length(what)/n1)
    if (!is.na(filen)) {ncols=n2;nlines=n1
    if (grepl(".jp",filen)) jpeg(filen, width = 480*3, height = 480*3,res=300)
    if (grepl(".pdf",filen)) pdf(filen)
    } else {ncols=n1;nlines=n2}
    par(mfrow=c(nlines,ncols))
    if (is.numeric(whichCell)) {cells=whichCell} else if (is.character(whichCell)) {
      if (whichCell=="all") {cells=cellNumA(validationdata$rasK)}
      if (whichCell=="cells_sampled") {cells<-as.numeric(levels(as.factor(validationdata$tip_genotypes$Cell_numbers)))}
      if (grepl("KmoreThan",whichCell)) {cells=which(t(validationdata$rasK)>strsplit(whichCell,"Than")[[1]][2])}
      if (whichCell=="KmoreThan1") {cells=which(t(validationdata$rasK)>1)}
    } else {cells=cellNumA(validationdata$rasK)}
    title =c(pop_size="Populations sizes",
             sample_size="Samples sizes",
             transition = "transition",
             Fsimul ="F simul",
             FlinearSimul="F/(1-F) Simul",
             geneticDistSimul="Simulations",
             linFstepping="F/(1-F) Step",
             genetDistStepping="g dist step",
             genetDistIsland="Island",
             hitting_time = "Hitting time",
             FlinearDigraph = "F/(1-F) Digraph",
             FlinearUndigraph="F/(1-F) Undigraph",
             genetDistUndigraph="Undigraph",
             genetiDistDigraph="Digraph",
             graphCoalTimes="coalescence time",
             graphCoalModelMutDist="Markov chain",
             genetDistGraphStepping="Stepping Stone",
             genetDistAbsorbingMC="Absorbant MC")
    matrixNames <- c(pop_size="popSize",
                     sample_size="sampleSize",
                     transition = "transitionmatrice",
                     Fsimul ="Fst_Simul",
                     FlinearSimul="linearizedFstSimul",
                     geneticDistSimul="geneticDistSimul",
                     linFstepping="linearizedFst_theory_stepping_stone_model",
                     genetDistStepping="genetDistStepping",
                     hitting_time = "hitting_time",
                     FlinearDigraph = "linearizedFst_theory_digraph_model",
                     FlinearUndigraph="linearizedFst_theory_undigraph_model",
                     genetDistUndigraph="genetDistUndigraph",
                     genetiDistDigraph="genetDistDigraph",
                     graphCoalTimes="coalescence_times_mean",
                     graphCoalModelMutDist="coalModelMutDist",
                     genetDistGraphStepping="genetDistGraphStepping",
                     genetDistIsland="genetDistIsland",
                     genetDistAbsorbingMC="genetDist_absorbingMethod")
    Cells1 <-  list(pop_size=1:dim(validationdata$rasK)[1],
                    sample_size=1:dim(validationdata$rasK)[1],
                    transition = 1:dim(validationdata$transitionmatrice)[1],
                    Fsimul =cells,
                    FlinearSimul=cells,
                    geneticDistSimul=1:dim(validationdata$geneticDistSimul)[1],
                    linFstepping=cells,
                    genetDistStepping=cells,
                    hitting_time = cells,
                    FlinearDigraph = cells,
                    FlinearUndigraph=cells,
                    genetDistUndigraph=cells,
                    genetiDistDigraph=cells,
                    graphCoalTimes=cells,
                    graphCoalModelMutDist=cells,
                    genetDistGraphStepping=cells,
                    genetDistAbsorbingMC=cells,
                    genetDistIsland=cells)
    Cells2=Cells1
    Cells2[["pop_size"]] <- 1:dim(validationdata$rasK)[2]
    Cells2[["sample_size"]] <- 1:dim(validationdata$rasK)[2]
    for (mat in what){
      plot(raster(validationdata[[matrixNames[mat]]][Cells1[[mat]],Cells2[[mat]]]),main=title[mat])
    }
    if(!is.na(filen)) dev.off()
  }
  
  
  
)


setMethod(
  f="stepwise",
  signature="",
  definition=function(coaltable,initial_locus_value,stepvalue,locusnames)
  {
    coaltable$dir = 2*rbinom(n=prod(dim(as.table(coaltable[,"mut"]))),unlist(coaltable[,"mut"]),.5)-coaltable[,"mut"]
    coaltable$genotype <- is.na(coaltable$mut)
    coaltable[dim(coaltable)[1],"genotype"] <- initial_locus_value
    for(branch in rev(rownames(coaltable)[-dim(coaltable)[1]]))
    {
      coaltable[branch,"genotype"] <- as.matrix(coaltable[branch,"dir"]*stepvalue + coaltable[which(coaltable$coalescing==coaltable[branch,"new_node"]),"genotype"])
    }
    coaltable
  }
  #
  #
  #
  
)


setMethod(
  f="matrixes_forward",
  signature="",
  definition=function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, shapeDisp, pDisp, a_value_obs, a_value_att, file=NULL, mutationRate,nbLocus, initial_locus_value,indpercell)
  {
    nblayers =dim(donneesEnvironmentObs)[3]
    nCell = ncellA(donneesEnvironmentObs)
    Cell_numbers <- cellNumA(donneesEnvironmentObs)
    geneticObs = simulationGenet(donneesEnvironmentObs,pK,pR,shapesK,shapesR,mutationRate,nbLocus,initial_locus_value,shapeDisp,pDisp,nb_generations=5000,indpercel)
    finalGenetData = geneticObs[[1]]
    a_value_simul = geneticObs[[2]]
    migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
    transitionmatrice = transitionMatrixForward(r, K, migration= migrationM)
    land_size <- raster(matrix(K[1,]))
    extent(land_size) <- extent(donneesEnvironmentObs)
    list(land_size,transitionmatrice,a_value_simul,a_value_att)
  }
  
  
  
)


setMethod(
  f="Qbetween_2",
  signature="",
  definition=function(gDat)
  {
    gDatG <- gDat[,grep("ocus",colnames(gDat))]
    gDat3D <- OneCol23Dims(gDatG)
    (id_mat_prob(gDat3D[,,1],gDat3D[,,1])
      +id_mat_prob(gDat3D[,,2],gDat3D[,,2])
      +id_mat_prob(gDat3D[,,1],gDat3D[,,2])
      +id_mat_prob(gDat3D[,,2],gDat3D[,,1]))/4
  }
  
)

