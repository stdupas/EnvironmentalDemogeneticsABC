library(doParallel)
registerDoParallel(cores=2)

ParallelGibbs <- function(n=8) {



    startF = c(2, 12, 6, 18, 2, 8, 6, 12)
    startC = c(1, 9, 3, 21, 1, 9, 3, 9)
    start = rbind(startF, startC)
    indiceF = 100 
    indiceC = 100
    indice = rbind(indiceF, indiceC)
    thining = 2
    nbPar = length(start[1,])

    ndvC = NULL
    ndvF = NULL
    ndpC = NULL
    ndpF = NULL


    scaleF = c(2,2,2,2,2,2,2,2)
    scaleC = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)
    scale = rbind(scaleF, scaleC)
    
    

    #res = rep(list(rep(0,nbPar),0),2)
    for(i in 1:n) {

        print("0")

        res = foreach(chaine = 1:2, .combine=c) %dopar%{
            #tmp = (chaine-1)*2+1
            #res[tmp:(tmp+1)] = 
            oneChainGibbs(start[chaine,], scale[chaine,], nbPar, indice[chaine,], thining)
        }

        print("1")

        postF = res[[2]]
        postC = res[[4]]
        paramF = res[[1]]
        paramC = res[[3]]

        print("2")

        if(postF > postC) {
            startC = paramF
            startF = paramF
        } else {
            startC = paramC
            startF = paramF
        }

        print("3")

        ndvC = cbind(ndvC,startC)
        ndvF = cbind(ndvF,startF)
        ndpC = cbind(ndpC,postC)
        ndpF = cbind(ndpF,postF)

        cat(i,"\n")
    }

    return(list(ndvC,ndpC,ndvF,ndpF))



}

oneChainGibbs <- function(start, scale, nbPar, indice, thining) {

    start0 = start
    post0 = logPostDens(start0)
    
    maxProb = post0
    maxParam = start

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
            if(start1[1]>=start1[2] || start1[1]>=start1[3] || start1[5]>=start1[6] || start1[5]>=start1[7]) {
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
    }
    return(list(maxParam,maxProb))
}