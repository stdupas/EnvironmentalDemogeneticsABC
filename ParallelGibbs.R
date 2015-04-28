library(doParallel)
registerDoParallel(cores=2)

ParallelGibbs <- function(n=5, nbPar=12, files=FALSE) {
    
    ##########################
    #
    # PARAMETRES DE DEPART
    #
    # Si files est TRUE, on reprend depuis les derniers parametres
    if(files == TRUE) {
        startF = read.table(file="PARAMFROID.txt",header=TRUE)
        startC = read.table(file="PARAMCHAUD.txt",header=TRUE)
        nbLines = dim(startF)[1]
        startF = startF[nbLines,]
        startC = startC[nbLines,]
        
        # Sinon, on part de nouveaux parametres
    } else {
        startF = c(0.5, 10, 4, 20, 0.5, 10, 4, 10, 270, 320, 295, 1)
        startC = c(0.5, 10, 4, 20, 0.5, 10, 4, 10, 270, 320, 295, 1)
        header=c("K.pr.Xmin", "K.pr.Xmax", "K.pr.Xopt", "K.pr.Yopt",
                 "R.pr.Xmin", "R.pr.Xmax", "R.pr.Xopt", "R.pr.Yopt",
                 "R.tas.Xmin", "R.tas.Xmax", "R.tas.Xopt", "R.tas.Yopt")
        write(header, file="PARAMFROID.txt", ncolumns=nbPar, append=FALSE)
        write(header, file="PARAMCHAUD.txt", ncolumns=nbPar, append=FALSE)
        
        write(header, file="PARAMFROID_oneChain.txt", ncolumns=nbPar, append=FALSE)
        write(header, file="PARAMCHAUD_oneChain.txt", ncolumns=nbPar, append=FALSE)
    }
    start = rbind(startF, startC)
    
    ##########################
    #
    # NOMBRE ITERATION POUR CHAQUE CHAINE
    #
    indiceF = 250 
    indiceC = 250
    indice = rbind(indiceF, indiceC)
    
    ##########################
    #
    # ECHELLE DE CHAQUE CHAINE
    #
    scaleF = c(1,1,1,1,1,1,1,1,1,1,1,1)
    scaleC = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)
    scale = rbind(scaleF, scaleC)
    
    ##########################
    #
    # THINING
    #
    thining = 1
    
    ##########################
    #
    # DECLARATION DES VARIABLES A STOCKER
    #
    ndvC = NULL
    ndvF = NULL
    ndpC = NULL
    ndpF = NULL
    
    ##########################
    #
    # DEBUT DU GIBBS SAMPLING 
    #
    for(i in 1:n) {
        
        cat(i, ": debut\n")
        
        res = foreach(chaine = 1:2, .combine=c) %dopar%{
            oneChainGibbs(start[chaine,], scale[chaine,], nbPar, indice[chaine,], thining, chaine)
        }
        
        postF = res[[2]]
        postC = res[[4]]
        paramF = res[[1]]
        paramC = res[[3]]
        
        if(postF > postC) {
            startC = paramF
            startF = paramF
            cat("Echange!\n")
        } else {
            startC = paramC
            startF = paramF
        }
        
        write(startF, file="PARAMFROID.txt",ncolumns=nbPar, append=TRUE)
        write(startC, file="PARAMCHAUD.txt",ncolumns=nbPar, append=TRUE)
        
        start = rbind(startF, startC)
        
        ndvC = cbind(ndvC,startC)
        ndvF = cbind(ndvF,startF)
        ndpC = cbind(ndpC,postC)
        ndpF = cbind(ndpF,postF)
        
        cat(i,": fin\n")
    }
    return(list(ndvC,ndpC,ndvF,ndpF))
}

oneChainGibbs <- function(start, scale, nbPar, indice, thining, chaine) {
    
    start0 = start
    post0 = logPostDens(start0)
    
    maxProb = post0
    maxParam = start
    
    for(i in 1:indice){
        for(j in 1:nbPar){ 
            start1 = start0
            # On pioche une valeur de pas pour faire bouger les hyperparametres a partir de start0
            start1[j] = start0[j] + rnorm(1) * scale[j]
            
            ##################### ConquadraticSkewed1
            ##
            ##
            if(start1[1]>=start1[2] || start1[1]>=start1[3] || start1[5]>=start1[6] || start1[5]>=start1[7] || start1[9]>=start[10] || start[9]>=start[11]) {
                start1[j] = start0[j]
            }

            ##
            ##
            #####################
            
            
            # On calcule les posteriors avec les nouveaux hyperparametres
            post1 = logPostDens(start1)
            # On decide si on garde ou non les nouvelles valeurs
            # Les valeurs sont gardees si la valeur tiree aleatoirement dans la loi uniforme est plus petite que
            # la difference entre le posterior au temps i et le posterior au temps i-1
            # t est egale a 1 si les valeurs sont gardees, sinon 0
            t = runif(1) < min(1,exp(post1 - post0))
            # Si t = 1, on garde les nouvelles valeurs, sinon on garde les anciennes valeurs
            start0[j] = start1[j] *(t==1) + start0[j] *(t==0)
            post0 = post1 * (t==1) + post0 * (t == 0)          
            
            if((i%%thining == 0) && (maxProb < post0)) {
                maxParam = start0
                maxProb = post0
            }
        }
        if(chaine==1) {
            write(start0, file="PARAMFROID_oneChain.txt",ncolumns=8, append=TRUE)
        }
        else {
            write(start0, file="PARAMCHAUD_oneChain.txt",ncolumns=8, append=TRUE)
        }
        
    }
    return(list(maxParam,maxProb))
}