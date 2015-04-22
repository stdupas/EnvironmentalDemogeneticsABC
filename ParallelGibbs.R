ParallelGibbs <- function(n=20) {
    startF = c(2, 33, 9, 2, 33, 2)
    startC = c(2, 33, 9, 2, 33, 2)
    start = rbind(startF, startC)
    indiceF = 400
    indiceC = 400
    indice = rbind(indiceF, indiceC)
    thining = 2
    nbPar = length(start)

    scaleF = c(1,1,1,1,1,1)
    scaleC = c(0.2,0.2,0.2,0.2,0.2,0.2)
    scale = rbind(scaleF, scaleC)
    
    system.time(
        foreach(i=1:2) %dopar%{
            for(j in 1:3){
                f(c[j,i])
            }
            
        })
    
    res = rep(list(rep(0,nbPar),0),2)
    for(i in 1:n) {
        foreach(chaine = 1:2) %dopar%{
            tmp = (chaine-1)*2+1
            res[[tmp:(tmp+1)]] = oneChainGibbs(start[chaine,], scale[chaine,], nbPar, indice[chaine,], thining)
        }


        postF = res[[2]]
        postC = res[[4]]
        paramF = res[[1]]
        paramC = res[[3]]

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

    start0 = start
    post0 = logPostDens(start0)
    
    maxProb = post0
    maxParam = NULL

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
            
            if((i%%thining == 0) && (maxProb < post0)) {
                maxParam = start0
                maxProb = post0
            }
        }
    }
    return(list(maxParam,maxProb))
}