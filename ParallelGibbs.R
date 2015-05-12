library(doParallel)
registerDoParallel(cores=2)

ParallelGibbs <- function(n=5, nbPar=12, files=FALSE) {
    
    ##########################
    #
    # DECLARATION DES VARIABLES A STOCKER
    #
    ndvC = NULL
    ndvF = NULL
    ndpC = NULL
    ndpF = NULL

    allStartC = NULL
    allPostC = NULL
    allStartF = NULL
    allPostF = NULL

    ##########################
    #
    # PARAMETRES DE DEPART
    #
    # Si files est TRUE, on reprend depuis les derniers parametres
    if(files == TRUE) {
        readF = as.numeric(read.table(file="PARAMFROID.txt",header=TRUE))
        readC = as.numeric(read.table(file="PARAMCHAUD.txt",header=TRUE))
        
        startF = readF[,(1:nbPar)] 
        startC = readC[,(1:nbPar)] 
        postF = readF[,(nbPar+1)]
        postC = readC[,(nbPar+1)]
        
        ndvC = c(ndvC,startC)
        ndvF = c(ndvF,startF)
        ndpC = c(ndpC,postC)
        ndpF = c(ndpF,postF)

        nbLines = dim(readF)[1]
        startF = startF[nbLines,]
        startC = startC[nbLines,]
        recovery2 = readRDS("PARAM_recovery.RData")

        allStartC = rbind(allStartC, read.table("ALLCHAUD_param.txt", header=TRUE, sep=" "))
        allPostC = rbind(allPostC, read.table("ALLCHAUD_post.txt", header=TRUE, sep=" "))
        allStartF = rbind(allStartF, read.table("ALLFROID_param.txt", header=TRUE, sep=" "))
        allPostF = rbind(allPostF, read.table("ALLFROID_post.txt", header=TRUE, sep=" "))
        
        # Sinon, on part de nouveaux parametres
    } else {
        startF = c(2, 15, 8, 25, 2, 15, 8, 15, 290, 310, 300, 2)
        startC = c(2, 15, 8, 25, 2, 15, 8, 15, 290, 310, 300, 2)
        # startF = c(0.5, 10, 4, 20, 0.5, 10, 4, 10, 285, 305, 295, 1)
        # startC = c(0.5, 10, 4, 20, 0.5, 10, 4, 10, 285, 305, 295, 1)
        header=c("K.pr.Xmin", "K.pr.Xmax", "K.pr.Xopt", "K.pr.Yopt",
                 "R.pr.Xmin", "R.pr.Xmax", "R.pr.Xopt", "R.pr.Yopt",
                 "R.tas.Xmin", "R.tas.Xmax", "R.tas.Xopt", "R.tas.Yopt", "posteriors")
        write(header, file="PARAMFROID.txt", ncolumns=nbPar+1, append=FALSE)
        write(header, file="PARAMCHAUD.txt", ncolumns=nbPar+1, append=FALSE)
        saveRDS(recovery2, "PARAM_recovery.RData")
        
        header=c("K.pr.Xmin", "K.pr.Xmax", "K.pr.Xopt", "K.pr.Yopt",
                 "R.pr.Xmin", "R.pr.Xmax", "R.pr.Xopt", "R.pr.Yopt",
                 "R.tas.Xmin", "R.tas.Xmax", "R.tas.Xopt", "R.tas.Yopt")
        write.table(allStartF, file="ALLFROID_param.txt", append=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")
        write.table(allStartC, file="ALLCHAUD_param.txt", append=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")
        
        write.table(allPostF, file="ALLFROID_post.txt", append=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")
        write.table(allPostC, file="ALLCHAUD_post.txt", append=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")

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
    scaleF = c(1,2,2,2,1,2,2,2,4,4,4,1)
    scaleC = c(0.2,0.4,0.4,0.4,0.2,0.4,0.4,0.4,0.8,0.8,0.8,0.2)
    scale = rbind(scaleF, scaleC)
    
    ##########################
    #
    # THINING
    #
    thining = 1
    
    ##########################
    #
    # DEBUT DU GIBBS SAMPLING 
    #
    for(i in 1:n) {
        
        cat(i, ": debut\n")
        
        res = foreach(chaine = 1:2, .combine=c) %dopar%{
            oneChainGibbs(start[chaine,], scale[chaine,], nbPar, indice[chaine,], thining)
        }
        
        postF = res[[2]]
        postC = res[[6]]
        paramF = res[[1]]
        paramC = res[[5]]

        allSC = res[[7]]
        allPC = res[[8]]
        allSF = res[[3]]
        allPF = res[[4]]
        
        if(postF > postC) {
            startC = paramF
            startF = paramF
            cat("Echange!\n")
        } else {
            startC = paramC
            startF = paramF
        }

        start = rbind(startF, startC)
        
        # Ecriture des fichiers a la suite
        write(c(startF,postF), file="PARAMFROID.txt", ncolumns=nbPar, append=TRUE)
        write(c(startC,postC), file="PARAMCHAUD.txt", ncolumns=nbPar, append=TRUE)

        write.table(allSF, file="ALLFROID_param.txt", append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ")
        write.table(allSC, file="ALLCHAUD_param.txt", append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ")
        
        write.table(allPF, file="ALLFROID_post.txt", append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ")
        write.table(allPC, file="ALLCHAUD_post.txt", append=TRUE, col.names=FALSE, row.names=FALSE, sep=" ")

        # Recuperation des start et post pour le max
        ndvC = cbind(ndvC,startC)
        ndvF = cbind(ndvF,startF)
        ndpC = cbind(ndpC,postC)
        ndpF = cbind(ndpF,postF)

        # Recuperation de toutes les valeurs de parametres et posteriors
        allStartC = rbind(allStartC,allSC)
        allPostC = rbind(allPostC,allPC)
        allStartF = rbind(allStartF,allSF)
        allPostF = rbind(allPostF,allPF)
        
        cat(i,": fin\n")
    }
    return(list(ndvC,ndpC,ndvF,ndpF))
}



################################
##
##
##          TEMPORAIRE
##
##
## Commandes pour lancer seulement le oneChainGibbs
## 
# start = c(2, 15, 8, 25, 2, 15, 8, 15, 290, 310, 300, 2)
# scale = c(1,2,2,2,1,2,2,2,4,4,4,1)
# nbPar=12
# indice = 100
# thining = 1
##
##
##
##
##
###############################


oneChainGibbs <- function(start, scale, nbPar, indice, thining) {
    # Fonction faisant tourner un algorithme de Gibbs Sampling
    # Variables: 
    #           start: vecteur contenant les valeurs de depart des parametres
    #           scale: vecteur contenant une valeur d'echelle pour le pas de chaque parametre
    #           indice: nombre d'iterations de l'algorithme
    #           nPar: nombre de parametres a evaluer
    #           ndv: tableau (nbIterations x nbParametres) contenant la valeur des parametres a chaque iteration
    #           ndp: tableau (nbIterations x nbParametres) contenant la valeur des posteriors a chaque iteration
    #           post0: posteriors a l'iteration (i-1), logPostDens
    #           start0: valeurs des hyperparametres a l'iteration (i-1)
    #           start1: valeurs des hyperparametres a l'iteration (i)
    #           post1: posteriors a l'iteration (i) 

    start0 = start
    post0 = logPostDens(start0)
    
    maxProb = post0
    maxParam = start

    allStart = array(0, dim = c(indice, nbPar))
    allPost = array(0, dim = c(indice, nbPar))
    
    for(i in 1:indice){
        cat("\n", i, ":")
        for(j in 1:nbPar){
            cat("*") 
            start1 = start0
            # On pioche une valeur de pas pour faire bouger les hyperparametres a partir de start0
            start1[j] = start0[j] + rnorm(1) * scale[j]
            
            ##################### ConquadraticSkewed1
            ##
            ##
            if(start1[1]>=start1[2] || start1[1]>=start1[3] || start1[5]>=start1[6] || start1[5]>=start1[7] || start1[9]>=start1[10] || start1[9]>=start1[11]) {
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

            #####################
            ##
            ## On garde toutes les valeurs de post et start
            allStart[i,j] = start0[j]
            allPost[i,j] = post0         
            
            # On garde les valeurs avec le meilleur posterior
            if((i%%thining == 0) && (maxProb < post0)) {
                maxParam = start0
                maxProb = post0
            }
        }  
    }
    return(list(maxParam,maxProb,allStart,allPost))
}