rm(list=ls())
wd="/media/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # portable
wd="/media/dupas/1To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # fixe
wd="/home/legs/GraphPOP/" # portable
wd="/home/dupas/GraphPOP/" # fixe
wd="/home/arno/These/GraphPOP" # portable arno
setwd(wd)
source("graphPOP_0.114.R")

########### Parameters initialisation  ########### >>>>>>

###### Environmental data of temperature and precipitations

# Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900)) 
Data2 <- data.frame(BIO1=c(300,120,120,400),BIO12=c(2000,350,350,2900)) 
# Make raster stack with two layers according to the environmental variables of the dataframe
rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1),"BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)))

###### Genetic parameters :
N=1.5
mutation_rate=1E-4
initial_genetic_value=200

##### Niche Function :

# Concave quadratic skewed distribution parameters
# Shapes of the reaction norms for demographic variables (r and K) :
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
shapesr=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
# Parameters of the reaction norm are given by pK and pr
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
pr = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))

# Linear model + constant model parameters
shapesK=c(BIO1="linearPositive",BIO12="linearPositive")
shapesr=c(BIO1="constant",BIO12="constant")
pK = matrix(c(100,0.01,300,0.001),nrow=2,ncol=2,dimnames=list(c("X0","slope"),c("BIO1","BIO12")))
pr = matrix(c(N,N),nrow=1,ncol=2,dimnames=list(c("Y"),c("BIO1","BIO12")))

##### Dispersion Function :

# Shape of the dispersion function :
shapeDisp="fat_tail1"
# Dispersion parameters of the dispersion function
pDisp = c(alpha=1/19,beta=1)

######################### end of parameters initialisation <<<<<<


######################### Coalescence Simulation Test >>>>>>>

##### Get the carrying capacity map :
rasK=rasterStack
values(rasK)= as.matrix(ReactNorm(X=values(rasterStack),p=pK,shapes=shapesK)[,"Y"])

###### Create genetic data 
# the genetic data is meaning less (only one genotype),
# but genetic data are not used. They will be modified according to the simulated coalescent
# we create as many individuals (lines) as K for each cell of the map 
# then we sample a few individuals for the coalescent
geneticData = CreateGenetArray(rasK=rasK, nb_locus=20, initial_locus_value=200, Option="full_1col_diploid",nind=4)          
geneticData = CreateGenetArray(rasK=rasK, nb_locus=20, initial_locus_value=sample((80:120)*2,6*20,replace=TRUE), Option="sample_1col_diploid", nind=3)
dim(geneticData)

##### Get the migration Matrix
migrationMatrix(rasterStack=rasK, shapeDisp=shapeDisp, pDisp=pDisp)

##### create a data and model lists to run simulations

data = list(geneticData=geneticData,
            rasterStack=rasterStack)

model = list(pK=pK, pr=pr,
             shapesK=shapesK, shapesr=shapesr,
             shapeDisp=shapeDisp, pDisp=pDisp,
             mutation_rate=1E-1, 
             initial_genetic_value=initial_genetic_value, 
             mutation_model="bigeometric",stepvalue=2,
             mut_param=c(p=.5,sigma2=4))

##### Simulate the coalescent
system.time(coalescent_simulated <- simul_coalescent(data,model))
coalescent_simulated
phylog_tree <- coalescent_2_phylog(coalescent=coalescent_simulated$coalescent)

##### Plot the coalescent
plot_coalescent(coalescent=coalescent_simulated$coalescent,genetic_table=coalescent_simulated$genetic_values,with_landscape=TRUE, rasK=rasK, legend_right_move=-.3)
plot_coalescent(coalescent=coalescent_simulated$coalescent,genetic_table=coalescent_simulated$genetic_values,rasK=rasK,legend_right_move=-.5)

##### New instance of class reaction norm


##### Appending summary stats reference table
new_reference_table(geneticData,Distance="Goldstein")
variables=c("pK.Xmin","pK.Xmax","pK.Xopt","pK.Yopt")
fill_reference_table <- function(geneticData=geneticData,Distance=Distance,
                                 rasterStack=rasterStack,
                                 pK=pK, pr=pr,
                                 shapesK=shapesK, shapesr=shapesr,
                                 shapeDisp=shapeDisp, pDisp=pDisp,
                                 mutation_rate=1E-1, 
                                 initial_genetic_value=initial_genetic_value, 
                                 mutation_model="tpm",stepvalue=2,
                                 mut_param=c(p=.5,sigma2=4))
  
######################### end of Coalescence Simulation Test <<<<<<<<

 
####### PSEUDO TRASH >>>>>>>>>>

###
setwd(wd)
source("graphPOP_0.113.R")
#wd="/media/dupas/LACIE SHARE/IRD/ETUDIANTS/2014/Vasundra/R/" # fixe
Data2 <- na.omit(read.table("Data.csv")[c("BIO1","BIO12")])
Data2 <- data.frame(BIO1=rep((11:40)*10,10),BIO12=rep((0:9)*400,each=30))
rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=30,ncol=10),xmn=0,xmx=10,ymn=0,ymx=30),"BIO12"=raster(matrix(Data2$BIO12,nrow=30,ncol=10),xmn=0,xmx=10,ymn=0,ymx=30)))
Data2 <- values(rasterStack)
pK = matrix(c(100,500,300,0,10,10,300,3000,2500,0,10,10),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
pr = matrix(c(100,500,300,0,10,10,300,3000,2500,0,10,10),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
shapesr=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
pDisp = c(3.6,1.6)
K = ReactNorm(values(rasterStack),pK,shapesK)[,"Y"]
r = ReactNorm(values(rasterStack),pr,shapesr)[,"Y"]
Data2 <- cbind(Data2,K)
rasK <- rasterStack;values(rasK) <- subset(Data2,select=K)
method="fat_tail1"
migrationM <- migrationMatrix(rasterStack,method, pDisp)
transitionmatrice = transitionMatrixBackward(r, K, migration= migrationM)
plot(raster(transitionmatrice))
geneticData = CreateGenetArray(rasK, 20,200,Option="full_pop")
geneticData = geneticData[sample(dim(geneticData)[1],100),];row.names(geneticData)=1:100
dim(geneticData)
coalescent_simulated <- simul_coalescent(geneticData,transitionmatrice,dimGeneticData,rasK)
coalescent_2_newick(coalescent_simulated)







transitionmatriceF =transitionMatrixForward(r, K, migration= migrationM)
donneesEnvironmentObs = stack(list(BIO1=raster(matrix(log(sample((5:50)*10,10)),nrow=5,ncol=2),xmn=0,xmx=2,ymn=0,ymx=5),
                                   BIO12=raster(matrix(log(sample((20:300)*10,10)),nrow=5,ncol=2),xmn=0,xmx=2,ymn=0,ymx=5)))

rasterStackReactNormInit <- rasterStack; values(rasterStackReactNormInit) <- cbind(values(rasterStack),r,K,d=.9,Init=c(rep(0,33),1,rep(0,266)))
plot(rasterStackReactNormInit)
dim( rasterStack)
dynamics <- rasterStack; values(dynamics) <- subset(values(rasterStackReactNormInit),select = Init)
plot(dynamics)
valdyn <- values(dynamics)[,rep(1,100)];dim(valdyn);colnames(valdyn)=paste(colnames(valdyn),1:100,sep="")
values(dynamics) <- valdyn

p = list(K=list(shape="linear0",priorshape = "unif", pK=matrix(c(50,300,300,600,10,10,100,250,250,600,10,10),nrow=6,ncol=2,dimnames=list(c("KXY0Min","KXY0Max","KXYmaxMin","KXYmaxMax","KYmaxMin","KYmaxMax"),c("BIO1","BIO12")))),
         r=list(shape="linear0",priorshape = "unif", pr=matrix(c(50,300,300,600,10,10,100,250,250,600,10,10),nrow=6,ncol=2,dimnames=list(c("rXY0Min","rXY0Max","rXmaxMin","rXmaxMax","rYXmaxMin","rYXmaxMax"),c("BIO1","BIO12"))),
                Tgen = list(shape="envelin",priorshape = "unif", pr=matrix(c(50,300,300,600,10,10),nrow=6,ncol=2,dimnames=list(c("rX00Min","rX00Max","rXmaxMin","rXmaxMax","rYXmaxMin","rYXmaxMax"),c("BIO1"))))))

matrix(c(100,500,300,0,10,10,100,500,300,0,10,10,100,500,60,30,300,3000,2500,0,10,10,300,3000,2500,0,10,10,300,3000,40,40),nrow=16,ncol=2,dimnames=list(c("KXmin","KXmax","KXopt","KYxmin","KYxmax","Yopt","rXmin","rXmax","rXopt","rYxmin","rYxmax","Yopt","GtimeXmin","GTimeXmax","GtimeYXmin","GTimeYXmax"),c("BIO1","BIO12")))
p_est = matrix(c(rep(TRUE,28),rep(FALSE,4)),nrow=16,ncol=2,dimnames=list(c("KXmin","KXmax","KXopt","KYxmin","KYxmax","Yopt","rXmin","rXmax","rXopt","rYxmin","rYxmax","Yopt","GtimeXmin","GTimeXmax","GtimeYXmin","GTimeYXmax"),c("BIO1","BIO12")))
pr = matrix(c(100,500,300,0,10,10,300,3000,2500,0,10,10),nrow=6,ncol=2,dimnames=list(c("rXmin","rXmax","rXopt","rYxmin","rYxmax","Yopt"),c("BIO1","BIO12")))
ptG = c(tgMean=30, tgSD=10) # paramètes temps de génération


values(dynamics) <- cbind(dynamics,gridRepnDispFunction(values(dynamics),r,K,d=.9,migrationM,overlapping=FALSE))

plotmig <- function(pDisp)
{
migrationM <- migrationMatrix(rasterStack,method, pDisp)
plot(raster(migrationM))
}

plotmig <- function(pDisp,max=1000)
{
  plot(0:100*max/100,fat_tail1(0:100*max/100,pDisp))
  c(sum(fat_tail1(0:100*max/100,pDisp)),sum(0:100*max/100*fat_tail1(0:100*max/100,pDisp)))
}

coalescent_simulated <- simul_coalescent(geneticData,transitionmatrice,dimGeneticData,Rast_donnees_taille)



rasterStack=exp(donneesEnvironmentObs)
Data2$z=combineReactNorms(values(rasterStack),p,c("quadraticskewed","quadraticskewed"));wireframe(z~BIO1*BIO12,data=Data2,scale=list(arrows=FALSE)) # requires library lattice
Data2 = values(exp(donneesEnvironmentObs))
enveloppe(Data2,p)
quadratic(Data2,p)
quadraticskewed(Data2,p)
rangelinear(Data2,p)
ReactNorm(Data2,p,"enveloppe")
ReactNorm(Data2,p,"quadratic")
ReactNorm(Data2,p,"rangelin")
ReactNorm(Data2,p,"rangeloglin")
ReactNorm(Data2,p,"quadraticskewed")

ReactNorm(Data2,p,"loG")
X=Data2
mutationRate = 10^(-2);                             methodDisp = "gaussian"; ind_per_cell=10
nbLocus=10;                                         initial_locus_value=200

donneesEnvironmentObs = raster(matrix(log(c(5,6)),nrow=1,ncol=2),xmn=0,xmx=2,ymn=0,ymx=1)
donneesEnvironmentObs = raster(matrix(log(c(2,2,2)),nrow=1,ncol=3),xmn=0,xmx=3,ymn=0,ymx=1)
donneesEnvironmentObs = raster(matrix(log(c(2,3,3,2)),nrow=2,ncol=2),xmn=0,xmx=2,ymn=0,ymx=2)
donneesEnvironmentObs = raster(matrix(log(sample(30,size=100,replace=TRUE)),nrow=10,ncol=10),xmn=0,xmx=10,ymn=0,ymx=10)
donneesEnvironmentObs = raster(matrix(log(100)),xmn=0,xmx=1,ymn=0,ymx=1)
donneesEnvironmentObs = raster(matrix(log(c(2,1)),nrow=1,ncol=2),xmn=0,xmx=2,ymn=0,ymx=1)
donneesEnvironmentObs = raster(matrix(log(10),nrow=5,ncol=2),xmn=0,xmx=2,ymn=0,ymx=5)
donneesEnvironmentObs = raster(matrix(log(sample(10,10)),nrow=5,ncol=2),xmn=0,xmx=2,ymn=0,ymx=5)
                                   
                                   

pDisp = res(donneesEnvironmentObs)[1]/qnorm(0.999,0,1);     nbpDisp =nbpDisp(methodDisp);
pDisp = res(donneesEnvironmentObs)[1]/qnorm(0.99,0,1);     nbpDisp =nbpDisp(methodDisp);
alpha = 0
beta =1
plot(donneesEnvironmentObs)
nblayers =dim(donneesEnvironmentObs)[3]
nCell = ncell(donneesEnvironmentObs)
  nblayers =dim(donneesEnvironmentObs)[3]
  nCell = ncell(donneesEnvironmentObs)
  Cell_numbers <- 1:nCell
  donnees_taille <-K_Function(donneesEnvironmentObs, alpha, beta) # recuperation du nb de parents par cellule = Npop
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,methodDisp,pDisp)
  transitionmatrice = transitionMatrixBackward(Npop = donnees_taille, migration= migrationM)
  donnees_taille <-K_Function(donneesEnvironmentObs, alpha, beta) # recuperation du nb de parents par cellule = Npop
  Rast_donnees_taille <- donneesEnvironmentObs ; values(Rast_donnees_taille) <- donnees_taille[1,]
  geneticData = CreateGenetArray(Rast_donnees_taille, nbLocus,initial_locus_value,Option="full_pop")
  geneticData[,"Cell_number_init"] <- geneticData[,"Cell_numbers"]
  dimGeneticData = dim(geneticData)
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,methodDisp,pDisp)
  transitionmatrice = transitionMatrixBackward(Npop = donnees_taille, migration= migrationM)
coalescent_simulated <- simul_coalescent(geneticData,transitionmatrice,dimGeneticData,Rast_donnees_taille)
plot_coalescent(coalescent_simulated)
coalescent <- coalescent_simulated
a=NULL;b=NULL
cells=sample(dim(transitionmatrice)[2],2,replace=FALSE);cells
for (i in 1:1000)
{
  a[i] <- simul_commute(cells,transitionmatrice)
  b[i] <- simul_coocur(cells,transitionmatrice)
}
mean(a);mean(b)*4;


df <-data.frame(time=append(a,b),meth = c(rep("commute",length(a)),rep("coocur",length(b))));anova(lm(time~meth,data=df))
tree2 <- read.tree("ex.tre")

tips = NULL
internals = NULL
nodes = NULL
for (i in 1:length(coalescent))
{
  nodes = append(nodes,coalescent[[i]]$coalescing,coalescent[[i]]$new_node)
  internals = append(internals,coalescent[[i]]$new_node)
}
nodes = as.numeric(levels(as.factor(c(nodes,internals))));nodes = nodes[order(nodes)]
tips = nodes[!((nodes)%in%(internals))]
tree <- "(a:2,(b:1,c:1):1);";cat(tree,file="ex.tre",sep="\n")
tree <- read.tree("ex.tre");tree$edge=tree$edge[-(1:4),]
ages = NULL;ages[tips]=0
for (i in 1:length(coalescent))
{
  ages[coalescent[[i]]$new_node] <- coalescent[[i]]$time
}
Edge = 1
for (i in 1:length(coalescent))
{
  for(coalescing in as.numeric(coalescent[[i]]$coalescing))#coalescing = coalescent[[i]]$coalescing[1]
  {
    tree$edge <- rbind(tree$edge,as.numeric(c(coalescent[[i]]$new_node,coalescing)))
    tree$edge.length[Edge] <- ages[coalescent[[i]]$new_node]-ages[coalescing]
    Edge=Edge+1
  }
}
tree$Nnode=length(internals)
tree$tip.label = as.character(tips)


tree$tip.label;tree$edge;tree$edge.length;tree$Nnode
tree2$tip.label;tree2$edge;tree2$edge.length;tree2$Nnode



sum(values(Rast_donnees_taille))
times=NA
for (i in 1:100)
{
  coalescent_simulated <- simul_coalescent(geneticData,transitionmatrice,dimGeneticData,Rast_donnees_taille)
  times[i]<- coalescent_simulated[[length(coalescent_simulated)]]$time
}
c(mean(times),var(times))


wd="/home/dupas/Bureau/IRD/ETUDIANTS/2014/Vasundra/R/"
donneesEnvironmentObs <- Aggregate_and_adjust_raster_to_data(raster(paste(wd,envdir,envfiles,sep="")),release=read.table(paste(wd,genetfile,sep="")), recovery=read.table(paste(wd,genetfile,sep="")), extend_band_size=0.1, aggregate_index=8)
donneesEnvirRandom <- donneesEnvironmentObs
plot(donneesEnvironmentObs)
values(donneesEnvirRandom) <- runif(ncell(donneesEnvironmentObs),min(na.omit(values(donneesEnvironmentObs))),max(na.omit(values(donneesEnvironmentObs))))
beta = 1/4
alpha= -min(na.omit(values(donneesEnvironmentObs)))/4
# donneesEnvirRandom; p = c(alpha,beta,80,1.5); methodDisp="fat_tail1"; nblayers=1; nbpDisp=2; file="graph_random.jpg"; mutationRate=1E-3; nbLocus=20; initial_locus_value=200;nb_generations=5000
result_graph_random <- validation(donneesEnvirRandom, p = c(alpha,beta,80,1.5), methodDisp="fat_tail1", nblayers=1, nbpDisp=2, file="graph_random.jpg", mutationRate=1E-3, nbLocus=20, initial_locus_value=200,nb_generations=5000)
# p = c(alpha,beta,0.1,1.5); methodDisp="contiguous"; nblayers=1; nbpDisp=1; file="graph_random_step.jpg"; mutationRate=1E-3; nbLocus=20; initial_locus_value=200;nb_generations=5000;indpercell=30
result_graph_random_stepping <- validation(donneesEnvirRandom, p = c(alpha,beta,0.1,1.5), methodDisp="contiguous", nblayers=1, nbpDisp=1, file="graph_random_step.jpg", mutationRate=1E-3, nbLocus=20, initial_locus_value=200,nb_generations=5000,indpercell=30)
result_island <- validation(raster(matrix(1,nrow=1,ncol=100),xmn=0,xmx=100,ymn=0,ymx=1), p = c(0,log(10),1), methodDisp="island", nblayers=1, nbpDisp=1, file="island.jpg", mutationRate=1E-3,nbLocus=20, initial_locus_value=200,nb_generations=5000,indpercell=1)
result_stepping <- validation(raster(matrix(1,nrow=1,ncol=100),xmn=0,xmx=100,ymn=0,ymx=1), p = c(0,log(10),1), methodDisp="contiguous", nblayers=1, nbpDisp=1, file="stepping.jpg", mutationRate=1E-3,nbLocus=20, initial_locus_value=200,nb_generations=5000)
result_small <- validation(raster(matrix(c(max(na.omit(values(donneesEnvironmentObs))),min(na.omit(values(donneesEnvironmentObs)+1))),nrow=1,ncol=5),xmn=0,xmx=5,ymn=0,ymx=1),  p = c(alpha,beta,80,1.5), methodDisp="contiguous", nblayers=1, nbpDisp=1, file="stepping.jpg", mutationRate=1E-3,nbLocus=20, initial_locus_value=200,nb_generations=5000)
# donneesEnvironmentObs=raster(matrix(c(max(na.omit(values(donneesEnvironmentObs))),min(na.omit(values(donneesEnvironmentObs)+1))),nrow=1,ncol=5),xmn=0,xmx=5,ymn=0,ymx=1);  p = c(alpha,beta,80,1.5); methodDisp="contiguous"; nblayers=1; nbpDisp=1; file="stepping.jpg"; mutationRate=1E-3;nbLocus=20; initial_locus_value=200;nb_generations=5000
donneesEnvironmentObs=raster(matrix(1,nrow=1,ncol=100),xmn=0,xmx=10,ymn=0,ymx=1); p = c(0,log(10),0.1); methodDisp="contiguous"; nblayers=1; nbpDisp=1; file="stepping.jpg"; mutationRate=1E-3;nbLocus=20; initial_locus_value=200;nb_generations=5000


# contiguous modele: dispersion de proche en proche populations de taille 10, taux de migration de 0.1, a-value=distance/(4Nm)=distance/4
#
#
donneesEnvironmentObs <- raster(matrix(1,nrow=1,ncol=100),xmn=0,xmx=100,ymn=0,ymx=1)
alpha=0;                          beta = log(10);
nbLocus = 10;                     initial_locus_value = 200;      mutationRate = 10^(-4);
methodDisp = "contiguous";          pDisp = 0.1;              nbpDisp =nbpDisp(methodDisp);
nblayers =dim(donneesEnvironmentObs)[3]
nCell = ncell(donneesEnvironmentObs)


geneticObs = simulationGenet(donneesEnvironmentObs,alpha, beta, mutationRate,nbLocus, initial_locus_value,methodDisp,pDisp,nb_generations=5000)
finalGenetData = geneticObs[[1]]
a_value_obs = geneticObs[[2]]
a_value_theory_stepping_stone_model = distanceMatrix(donneesEnvironmentObs)/(4*0.05)
plot_matrixes_forward(donneesEnvironmentObs, alpha, beta, methodDisp, pDisp, a_value_obs, a_value_theory_stepping_stone_model,file=paste(wd,"stepping stone model.jpg",sep=""))
plot_matrixes_forward(donneesEnvironmentObs, alpha, beta, methodDisp, pDisp, a_value_obs, a_value_theory_stepping_stone_model,file=NULL)
Cell_numbers <- 1:ncell(donneesEnvironmentObs)
a_value_theory_graph_model <- expect_a_value(donneesEnvironmentObs,p=c(alpha,beta,pDisp),Cell_numbers,nbpDisp=1,nblayers=1,methodDisp)


#a_value_obs[a_value_obs<0]<-0 # (pour plus de clarté, on met à 0 les a_value_obs <0, c'est juste pour comparer pour l'instant)
#Stepping stone one dimensionnal model : slatkin 
#coords_final= xyFromCell(donneesEnvironmentObs, c(1:nrow(finalGenetData)))
#distFinal = as.matrix(dist(coords_final))
#dist_genet = distFinal[finalGenetData[,"Cell_numbers"],finalGenetData[,"Cell_numbers"]]
#x = dist_genet/(4*100*0.1)


# Island modele: homogen dispersal to cells M=1; Fst=0.2
donneesEnvironmentObs <- raster(matrix(1,nrow=1,ncol=100),xmn=0,xmx=10,ymn=0,ymx=1)
alpha=0;                        beta = log(10);
nbLocus = 20;                   initial_locus_value = 200;        mutationRate = 10^(-4);
methodDisp = "island";          pDisp = 0.1;                 nbpDisp =nbpDisp(methodDisp);
nblayers =dim(donneesEnvironmentObs)[3]
nCell = ncell(donneesEnvironmentObs)
geneticObs = simulationGenet(donneesEnvironmentObs,alpha, beta, mutationRate,nbLocus, initial_locus_value,methodDisp,pDisp,nb_generations=5000)
finalGenetData = geneticObs[[1]]
a_value_obs = geneticObs[[2]]
a_value_att = matrix(0.2,nrow=100,ncol=100)-0.2*diag(100)
plot_matrixes_forward(donneesEnvironmentObs, alpha, beta, methodDisp, pDisp, a_value_obs, a_value_att,file=paste(wd,"island model.jpg",sep=""))
#a_value_obs[a_value_obs<0]<-0
#par(mfrow=c(2,2))
#plot(raster(a_value_obs))

#x = 1/4*nCell*pDisp

# 

# fichier format genepop
coords_final= xyFromCell(donneesEnvironmentObs, finalGenetData[,"Cell_numbers"])
genotype_final = finalGenetData[,grep("Locus", colnames(finalGenetData), fixed = T)]
gl2gp(coords_final,genotype_final,"fichier_genepop.txt") #geneland to genepop fichier
gpop <- readLines("fichier_genepop.txt")
gpop <- gpop[gpop!="Pop"]
starting <- min(grep("sample",gpop))
gpop2 <- gpop[1:(starting-1)]
for (i in grep("sample",gpop))#i=9
{
  gpop2[starting+2*(i-starting)] <- "POP"
  gpop2[starting+2*(i-starting)+1] <- gsub(")","",strsplit(gpop[i],"[(]")[[1]][2])
}
writeLines(gpop2,"fichier_genepop.txt")

# Modèle simple
mutationRate = 10^(-2);                             methodDisp = "gaussian"; ind_per_cell=10
nbLocus=10;                                         initial_locus_value=200
donneesEnvironmentObs = raster(matrix(log(c(5,5)),nrow=1,ncol=2),xmn=0,xmx=2,ymn=0,ymx=3)
donneesEnvironmentObs = raster(matrix(log(c(10,2,10)),nrow=3,ncol=1),xmn=0,xmx=1,ymn=0,ymx=3)
donneesEnvironmentObs = raster(matrix(log(c(10,2,10,10,2,10)),nrow=3,ncol=2),xmn=0,xmx=2,ymn=0,ymx=3)
donneesEnvironmentObs = raster(matrix(log(c(rep(50,5),rep(20,5),rep(10,5),rep(20,5),rep(50,5))),nrow=5,ncol=5),xmn=0,xmx=5,ymn=0,ymx=5)

pDisp = res(donneesEnvironmentObs)[1]/qnorm(0.99,0,1);     nbpDisp =nbpDisp(methodDisp);

alpha = 0
beta =1
plot(donneesEnvironmentObs)
nblayers =dim(donneesEnvironmentObs)[3]
nCell = ncell(donneesEnvironmentObs)
result_small <- validation(donneesEnvironmentObs,  p = c(alpha,beta,pDisp[1],pDisp[2]), methodDisp=methodDisp, nblayers=1, nbpDisp=1, file="gaussian_strip.jpg", mutationRate=mutationRate,nbLocus=20, initial_locus_value=200,nb_generations=5000)

geneticObs = simulationGenet(donneesEnvironmentObs,alpha=0, beta = 1,mutationRate,nbLocus, initial_locus_value,methodDisp,pDisp,nb_generations=5000)
finalGenetData = geneticObs[[1]]
a_value_obs = geneticObs[[2]]


## Modele environnemental avec separation par une "vallée":
mutationRate = 10^(-4);                             methodDisp = "gaussian";
nbLocus=10;                                         initial_locus_value=200
pDisp = res(donneesEnvironmentObs)[1]/qnorm(0.975,0,1);     nbpDisp =nbpDisp(methodDisp);
donneesEnvironmentObs = rasterCrop
alpha = 0
beta =1

x=log(c(4,4,1,1,1,5,5,4,4,1,1,1,5,5,2,1,1,1,1,1,2,2,1,1,1,1,1,2,2,2,1,1,1,3,3,2,2,1,1,1,3,3))
values(donneesEnvironmentObs) = x
plot(donneesEnvironmentObs)

nblayers =dim(donneesEnvironmentObs)[3]
nCell = ncell(donneesEnvironmentObs)
geneticObs = simulationGenet(donneesEnvironmentObs,alpha=0, beta = 1,mutationRate,nbLocus, initial_locus_value,methodDisp,pDisp,nb_generations=5000)
finalGenetData = geneticObs[[1]]
a_value_obs = geneticObs[[2]]



alpha=0; beta1=0; beta2=0; disp1=0; disp2=0
for(i in seq(0,0.5,0.01)){
  initial_alpha= c(i,1.5,0.5,0.6,86)
  initial_beta1 = c(0,i,0.5,0.6,86)
  initial_beta2 = c(0,1.5,i,0.6,86)
  initial_disp1 = c(0,1.5,0.5,i,86)
  initial_disp2 = c(0,1.5,0.5,0.6,i)
  alpha[i*100] =ssr(initial_alpha)
  beta1[i*100] =ssr(initial_beta1) 
  beta2[i*100] =ssr(initial_beta2) 
  disp1[i*100] =ssr(initial_disp1)
  disp2[i*100] =ssr(initial_disp2)
  print(i)
}

par(mfrow=c(3,2))

plot(alpha, type="l")
plot(beta1, type="l")
plot(beta2, type="l")
plot(disp1, type="l")
plot(disp2, type="l")
par(mfrow=c(1,1))

plot(alpha)
plot(beta1)
plot(beta2)
plot(disp1)
plot(disp2)
initial = c(0,1,0.02,2,1)
fct_erreur_calc = nlm(f = ssr , p = initial, hessian = FALSE)
ssr(initial)
optim(initial , fn = ssr, hessian=F )
ssr(initial)

envdir="environment/"
#envfiles=c("bio8_16.bil","tmean8_16.bil","prec8_16.bil", "alt_16.bil")
#envfiles=c("bio8_16.bil","tmean8_16.bil")
envfiles=c("bio8_16.bil")
genetfile="/donnees_ird/DataCamargue.csv"
aggregate_factor=16
extention_of_zone_of_interest=0.1
projectionsystem="+proj=longlat +datum=WGS84" #set projection system of your maps

# Load maps in a raster stack and assign the projection system to the rasterStack
rasterStack <- stack(paste(wd,envdir,envfiles,sep="")) #Stack raster objects
proj4string(rasterStack) <- CRS(projectionsystem)

# Load genetic data:
donnees_Camargue <- read.table(paste(wd,genetfile,sep=""), header=T)
donnees_Camargue = SpatialPointsDataFrame(donnees_Camargue[,c("longitude","latitude")],donnees_Camargue, proj4string = CRS(projectionsystem))
