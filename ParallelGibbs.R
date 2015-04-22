ParallelGibbs <- function(n=20) {
    startF = rbind(2, 33, 9, 2, 33, 2)
    startC = rbind(2, 33, 9, 2, 33, 2)
    indiceF = 400
    indiceC = 400
    thining = 2
    nbPar = length(start)

    scaleF = c(1,1,1,1,1,1)
    scaleC = c(0.2,0.2,0.2,0.2,0.2,0.2)

    for(i in 1:n) {
        cf = parallel(oneChainGibbs(startF, scaleF, nbPar, indiceF, thining), name="froid")
        cc = parallel(oneChainGibbs(startC, scaleC, nbPar, indiceC, thining), name="chaud")

        res = collect(list(cf,cc), wait=TRUE)

        postF = res$froid[1,2]
        postC = res$chaud[1,2]
        paramF = res$froid[,1]
        paramC = res$chaud[,1]

        if(postF > postC) {
            startC = paramF
            startF = paramF
        } else {
            startC = paramC
            startF = paramF
        }

        cat(i,"\n")
        cat("Chaud:",paramC,"\n")
        cat("Froid:",paramF,"\n")

    }

}

oneChainGibbs <- function(start, scale, nbPar, indice, thining) {
    ndv = array(0, dim=c(indice, nbPar))
    ndp = array(0, dim = c(indice, nbPar))

    start0 = start
    post0 = logPostDens(start0)
    
    indiceMax = 0
    postMax = post0

    for(i in 1:indice){
        cat("\n", i, ":")
        for(j in 1:nbPar){ 
            cat("*")
            start1 = start0
            # On pioche une valeur de pas pour faire bouger les hyperparametres a partir de start0
            start1[j] = start0[j] + rnorm(1) * scale[j]
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
            ndv[i,j] = start0[j]
            ndp[i,j] = post0
            
            if((i%%thining == 0) && (postMax < post0)) {
                indiceMax = i
                postMax = post0
            }
        }

        # if(i%%thining == 0) {
        #     cat("\n",start0,"\n")
        # }
    }
    return(list(ndv,ndp,indiceMax,postMax))
}