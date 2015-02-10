####
#### STATISTIC MODEL : INFERENCE OF ECOLOGICAL MODEL FROM GENETIC DATA
####

##########################################################################
############## Set your working directory and files to load ##############
##########################################################################

# matrice : colonne
#raster : ligne



####################################################
##### BACKWARD MODEL FUNCTIONS AND EXECUTIONS ######
#####  SIMULATION OF PREDICTED GENETIC DATA   ######
####################################################


degree2km = function(rasterStack){
  # Function to get spatial resolution in km from a rasterStack
  #
  # Args:
  #   rasterStack: the rasterStack from which to obtain the resolution
  #
  # Returns:
  #   The spatial resolution in km from the rasterStack 
  x_origin = ((xmin(rasterStack)+xmax(rasterStack))/2) #longitude origin
  y_origin = ((ymin(rasterStack)+ymax(rasterStack))/2) #latitude origin
  x_destination = (x_origin + xres(rasterStack)) #longitude of destination point
  y_destination = (y_origin + yres(rasterStack)) #latitude of destination point
  
  dist_degree <- acos(sin(x_origin)*sin(x_destination)+cos(x_origin)*cos(x_destination)*cos(y_origin-y_destination))
  dist_km = dist_degree * 111.32
  dist_km
}


Aggregate_and_adjust_raster_to_data <- function(Envir_raster_stack,release,recovery,extend_band_size,aggregate_index)
{
  # Change resolution and extent of environmental stacked layers according to data geographic range and extension zone outside geographic range of data
  #
  # Args:
  #   Envir_raster_stack: raster file 
  #   release: release points file (columns "X" and "Y" as longitude nd latitude)
  #   recovery: recovery points file (columns "X" and "Y" as longitude nd latitude)
  #
  # Returns:
  #   The transformed rasterStack
  samples <- SpatialPoints(rbind(na.omit(release[,c("X","Y")]),na.omit(recovery[,c("X","Y")])))
  if (aggregate_index > 1) {Envir_raster_stack <- aggregate(crop(Envir_raster_stack,extent(samples)+extend_band_size), fact=aggregate_index, fun=mean, expand=TRUE, na.rm=TRUE)} else {
    Envir_raster_stack <- crop(Envir_raster_stack,extent(samples)+extend_band_size)
  }
  Envir_raster_stack
}

Show_Niche <- function(BBox,nb_points,p,shapes=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")) # non terminé
{
  # Allow to visualize two dimensional niche function 
  #
  # Args:
  #   BBox: pounding box of variable values (two variables) data.frame columns as variable names, lines as c("Min","Max") 
  #   nb_points: number of points to draw between min and max for each variable
  #   p: parameter values of the reaction norm for each variable as column
  #   shapes: shapes of the reaction norms for each variable in a vector
  #
  # Returns:
  #
  # Example
  # BB = matrix(c(100,400,200,3200),nrow=2,ncol=2,dimnames=list(c("Min","Max"),c("BIO1","BIO12")))
  # p = matrix(c(100,500,300,0,10,10,300,3000,2500,0,20,20),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
  # Shapes=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
  # Show_Niche(BB,nb_points=c(12,18),p,shapes=Shapes)
  
  Data = as.data.frame(matrix(NA,nrow=1,ncol=length(shapes)));colnames(Data)=colnames(p);Data=Data[-1,]
  n=rep(1,length(shapes))
  Var1=NULL
  for(i in 1:nb_points[1])
  {
    Var1=append(Var1,rep(i,nb_points[2]))
  }
  Var2 = rep(1:nb_points[2],nb_points[1])
  Data = as.matrix(data.frame(Var1=Var1,Var2=Var2));colnames(Data)=colnames(p)
  Data = BB[rep("Min",dim(Data)[1]),]+(Data-1)*(BB[rep("Max",dim(Data)[1]),]-BB[rep("Min",dim(Data)[1]),])/matrix(nb_points-1,nrow=dim(Data)[1],ncol=dim(Data)[2],byrow=TRUE)
  rownames(Data)=1:dim(Data)[1];Data=as.data.frame(Data)
  form = as.formula(paste("z~",paste(names(shapes),collapse="*"),sep=""))
  Data[,"z"]=ReactNorm(Data,p,shapes)[,"Y"]
  wireframe(form,data=Data,scales=list(arrows=FALSE)) # requires library lattice
}

# populationSize: uses K_Function to obtain the populationsize landscape raster
#
#
populationSize <- function(donneesEnvironmentObs, p, shapes)
{
  # Give population size according to a landscape raster.
  #
  # Args:
  #   donneesEnvironmentObs:
  #   p:
  #   shapes:
  #
  # Returns:
  #   Population size
  populationSize <- donneesEnvironmentObs
  values(populationSize) <- ReactNorm(valules(donneesEnvironmentObs), p, shapes)[1,]
  populationSize
}
  
distanceMatrix <- function(rasterStack){
  # (optional) distanceMatrix return distance between all cells of raster
  # get x and y coordinates for each cell of raster object put in parameters
  coords = xyFromCell(rasterStack, 1:length(values(rasterStack[[1]])), spatial=FALSE)
  distance = as.matrix(dist(coords)) # distance matrix of coordinates
  return(distance)
}



#prior

# simulation forward of population sizes across time
forward_simul_landpopsize <- function(N0,p, migration)
{
  
}
# laplaceMatrix returns Laplacian matrix from transition matrix
laplaceMatrix <- function(transitionMatrix){
  matrixD = matrix(0,nrow = dim(transitionMatrix), ncol = dim(transitionMatrix))
  diag(matrixD) = 1 # diagonal equals to 1
  laplacianMatrix = matrixD - transitionMatrix
  laplacianMatrix[is.na(laplacianMatrix)]<-0 # replace NA by 0
  laplacianMatrix
  
}

# Calcul of resistance between two points of the graph
# with the Moore-Penrose generalized inverser matrix.
# ref : Bapat et all, A Simple Method for Computing Resistance Distance (2003)
# ref : Courrieu, Fast Computation of Moore-Penrose Inverse Matrices (2005)
resistDist <- function(laplacianMatrix){
  inverseMP = ginv(laplacianMatrix) # generalized inverse matrix  (Moore Penrose)
  diag = diag(inverseMP) # get diagonal of the inverse matrix
  mii = matrix(diag, nrow =dim(inverseMP), ncol = dim(inverseMP))
  mjj = t(mii)
  mij = inverseMP
  mji = t(mij)
  commute_time = mii + mjj - mij - mji
  commute_time
}


# Calcul of genetic distance from resistance 
geneticDist <- function(commute_time, popSize){
  #genetic_dist = commute_time / (8* popSize)
  genetic_dist = commute_time / (8* (sum(popSize)/(dim(popSize)[1]*dim(popSize)[2])))
  genetic_dist
}

# MAIN EXECUTIONS DES FONCTIONS PERMETTANT D'OBTENIR IN FINE LES DONNEES GENETIQUES PREDITES
#rasterCrop = Aggregate_and_adjust_raster_to_data(raster(paste(wd,envdir,envfiles,sep="")),release=read.table(paste(wd,genetfile,sep="")), recovery=read.table(paste(wd,genetfile,sep="")), extend_band_size=1, aggregate_index=aggregate_factor)
#plot(rasterCrop)
  

###################################################
##### FORWARD MODEL FUNCTIONS AND EXECUTIONS ######
#####  SIMULATION OF OBSERVED GENETIC DATA   ######
###################################################


#Plot genetic data in environmental data observed 
#Genetic data is turn into a Spatial Pixel Data Frame 
#Mettre des couleurs en fonction du nombre d'individu
plotGeneticData = function(geneticData, EnvironmentalDataObserved){
  colnames(geneticData)[1:2] <- c("x","y")
  geneticData = SpatialPixelsDataFrame(points = geneticData[,c("x","y")], data = geneticData[,])
  plot(EnvironmentalDataObserved[[1]])
  plot(geneticData, add = T)
}

# Function that computes population size distribution moments in a grid from one generation to the other
# N population sizes of the parent genration
# r growth rates
# K carrying capacities
# d death rates
# migration transition matrix between cells of the grid
# ptG : parameters of generaiton time model

gridRepnDispFunction <- function(dynamics,r,K,d=.9,ptG, migration,overlapping=TRUE)
{
  # values(dynamics)[,dim(values(dynamics))[2]] is value at previoud day
  # d is mortality
  Nt = values(dynamics)[,dim(values(dynamics))[2]]*(1-d) + r*N*(K-N/K)
  esperance[K==0] <- 0  
}

# Function that combine reproduction, dispersion and mutation for a given genetic data
repnDispMutFunction <- function(geneticData, dimGeneticData, mutationRate, transitionmatrice){
  # Calcul for reproduction and dispersion 
  # random choice of individuals and at the same time of their target cells in the transition matrix
  ncell_transition <- dimGeneticData[1]*dim(transitionmatrice)[1]
  transition_celnu <- sample(ncell_transition,dimGeneticData[1],replace=FALSE,transitionmatrice[geneticData[,"Cell_numbers"],])
  transition_col <- ceiling(transition_celnu/dimGeneticData[1])
  transition_line <- transition_celnu%%(dimGeneticData[1]);transition_line[transition_line==0]<-dimGeneticData[1]
  cell_numbers_sampled <- geneticData[transition_line,"Cell_numbers"]
  geneticData <- geneticData[transition_line,]
  geneticData[,"Cell_numbers"] <- transition_col
  locusCols = grep("Locus", colnames(geneticData))
  step = 2
  mu = mutationRate # mutation rate
  liability = runif(prod(dimGeneticData), 0, 1) # mutation liability
  liability = as.data.frame(matrix(liability, ncol = length(locusCols), nrow = dimGeneticData[1]))
  geneticData[,locusCols] = geneticData[,locusCols] + ((liability<mu/2)*step - (liability>(1-mu/2))*step) 
  #print(c("mutrat",(sum(liability<mu/2)+sum(liability>(1-mu/2)))/(length(grep("Locus", colnames(geneticData)))*dimGeneticData[1])))
  geneticData
}

# Function that combine reproduction, dispersion and mutation for a given genetic data
repnDispMutFunction <- function(geneticData, dimGeneticData, mutationRate, transitionmatrice){
  # Calcul for reproduction and dispersion 
  # loop for individuals
  locusCols = grep("Locus", colnames(geneticData))
  for (individual in 1:dimGeneticData[1])
  { # we choose where the parent come from in the individual cell line probabilities of the backward transition matrice
    mothercell = sample(nCell, 1,,transitionmatrice[geneticData[individual,"Cell_numbers"],])
    # we chose the parent among the individuals in this cell
    geneticline = sample(which(geneticData[,"Cell_numbers"]==mothercell),1)
    # we atribute the individual the genetic data of its mother
    geneticData[individual,locusCols] = geneticData[geneticline,locusCols]
  }
  step = 2
  mu = mutationRate # mutation rate
  liability = runif(prod(dimGeneticData), 0, 1) # mutation liability
  liability = as.data.frame(matrix(liability, ncol = length(locusCols), nrow = dimGeneticData[1]))
  geneticData[,locusCols] = geneticData[,locusCols] + ((liability<mu/2)*step - (liability>(1-mu/2))*step) 
  geneticData
}

#Function that calculate probability of identity of genes intra individual (at individuals level)
Qwithin_pair <- function(geneticData){
  matrix_pair = geneticData[,grep(".2", colnames(geneticData), fixed = T)]
  matrix_impair = geneticData[,grep(".1", colnames(geneticData), fixed = T)]
  Qw <- (matrix_pair == matrix_impair)
  Qw = rowMeans(Qw) # vector of probability of Qw for each individual
  Qw = matrix(Qw, ncol = dim(geneticData)[1], nrow = dim(geneticData)[1])
  Qw = (Qw+t(Qw))/2
  Qw
}

#Function that calculate probability of identity of genes intra individual (at population level)
Qwithin_pop <- function(geneticData){
  matrix_pair = geneticData[,grep(".2", colnames(geneticData), fixed = T)]
  matrix_impair = geneticData[,grep(".1", colnames(geneticData), fixed = T)]
  Qw <- (matrix_pair == matrix_impair)
  Qw = mean(Qw*1)
  Qw
}


# Fonction that calculates probability of identity of genes inter individual (between two individuals)
Qbetween <- function(geneticData, dimGeneticData){
  Qb = matrix(ncol = dimGeneticData[1], nrow = dimGeneticData[1]) #initialization of Qb as a matrix
  # A = genetic data with loci only
  A=as.matrix(geneticData[,grep("Locus",colnames(geneticData),fixed=T)])
  # On construit un tableau à plusieurs étages: il y a autant d'étage qu'il y a de locus ( les alleles passent en étage)
  A3 = aperm(array(A,dim=c(dim(A)[1],dim(A)[2],dim(A)[1])), c(1,3,2)) # permutation des colonnes et des etages
  B3 = aperm(A3, c(2,1,3)) # transposee de A3
  moy1 = colMeans(aperm(A3 == B3), dims = 1, na.rm = T)
  
  #Permutation des colonnes deux à deux des alleles de A pour calculer l'autre cas possible d'identite des alleles
  l= 1:dim(A)[2]
  Aprime= A[,c(matrix(c(l[2*floor(1:(length(l)/2))],l[2*floor(1:(length(l)/2))-1]), nrow= 2, ncol = length(l)/2, byrow = T))] # Permutation des colonnes
  # permutation et creation des etages pour Aprime ( on ne change qu'une seule des matrice (A3 / B3) et l'autre est inchangée pour comparer les alleles"complementaires"
  Aprime3 = aperm(array(Aprime,dim=c(dim(A)[1],dim(A)[2],dim(A)[1])), c(1,3,2))
  moy2 = colMeans(aperm(Aprime3 == B3), dims = 1, na.rm = T) # calcul moy dist pour les individus avec loci permutés Aprime3 et la transposée B3
  #Mean of distance between individuals: Qbetween
  Qb =(moy1 + moy2)/2
  Qb
}


# TEST: pour un nombre de generation donné, on teste la stabilité du a value
#Function test of stabilisation for a value
test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice){
## ref: Rousset et al. J Evol Biol,13  (2000) 58-62.
  vecteur_a_value <-c(0)
  for(i in 1: nb_generations){
    print(i)
    geneticData <- repnDispMutFunction(geneticData, dimGeneticData, mutationRate, transitionmatrice) 
    matrixQb = (1-Qwithin_pair(geneticData)+2*(Qwithin_pair(geneticData)-Qbetween(geneticData,dimGeneticData)))
    matrixQw = 2*(1-Qwithin_pop(geneticData))
    vecteur_a_value[i] = matrixQb/matrixQw-1/2
    vecteur_a_value[is.na(vecteur_a_value)] <-0
    if((i>90) && (i%%30 == 0)){
      if(var(vecteur_a_value[(i-30):i])> var(vecteur_a_value[(i-60):(i-30)]) 
         && var(vecteur_a_value[(i-30):i])> var(vecteur_a_value[(i-90):(i-60)])){
        return(list(geneticData, (matrixQb/matrixQw-1/2)))
        break
      }
    }
  }
}

test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice){
  ## ref: Rousset et al. J Evol Biol,13  (2000) 58-62.
  vecteur_a_value <-c(0)
  for(i in 1: nb_generations){
    print(i)
    geneticData <- repnDispMutFunction(geneticData, dimGeneticData, mutationRate, transitionmatrice) 
    Genotypes = geneticData[,grep("Locus", colnames(geneticData), fixed = T)]
    popmbrship=geneticData[,"Cell_numbers"]
    Fst = Fstat(Genotypes,nCell,popmbrship,ploidy=2)
    Fstlinear[i] <- Fst/(1-Fst)
    if((i>90) && (i%%30 == 0)){
      if(var(Fst[(i-30):i])> var(Fst[(i-60):(i-30)]) 
         && var(Fst[(i-30):i])> var(Fst[(i-90):(i-60)])){
        return(list(geneticData, Fst))
               break
      }
    }
  }
}



fstat = function(geneticData){
  Genotypes = geneticData[,grep("Locus", colnames(geneticData), fixed = T)]
  form <- as.formulae
  Pops = geneticData[,"Cell_numbers"]
  MeanPop = t(matrix((colSums(Genotypes)/dimGeneticData[1]),ncol=dimGeneticData[1],nrow=dimGeneticData[2]))
  VarInd = matrix(Genotypes^2 - MeanTot^2,ncol=dimGeneticData[2])
  VarInterPop = var(MeanPop)
  VarIntraPop = colSums(VarInd)/dimGeneticData[1]
  VarTot = VarInd
}

simul_commute <- function(cells=c(1,2),transitionmatrice)
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

simul_coocur <- function(cells=c(1,2),transitionmatrice)
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

combine.names = function(names1,names2)
{
  combination=NULL
  for (name1 in names1)
  {
    for (name2 in names2)
    {
      combination=append(combination,paste(name1,name2,sep="."))
    }
  }
combination
}


react_norm_param <- function(shapes)
{
  params = NULL
  for (shape in shapes)
  {
  params = append(params, switch(shape,
           constant = paste(names(rasterStack),".Y",sep=""),
           enveloppe = combine.names(names(rasterStack),c("Xmin","Xmax","Yopt")),
           envelin = combine.names(names(rasterStack),c("Yxmin","Yxmax","Xmin","Xmax")),
           envloglin = combine.names(names(rasterStack),c("Yxmin","Yxmax","Xmin","Xmax")),
           linear = combine.names(names(rasterStack),c("X0","slope")),
           linearPositive = combine.names(names(rasterStack),c("X0","slope")),
           conquadratic = combine.names(names(rasterStack),c("Xmin","Xmax","Yopt")),
           conquadraticsq = combine.names(names(rasterStack),c("Xmin","Xmax","Yopt")),
           conquadraticskewed = combine.names(names(rasterStack),c("Xmin","Xmax","Xopt","Yopt")),
           conquadraticskewedsq = combine.names(names(rasterStack),c("Xmin","Xmax","Xopt","Yopt"))
                                )
           )
  }
params
}

input_reaction_norm_shape_model <- function(demographic_parameter,names_envir_variables)
{
allshapes=  c("constant", "enveloppe", "envelin", "envloglin", "linear", "linearPositive", "conquadratic", "conquadraticsq", "conquadraticskewed", "conquadraticskewedsq")
shape=NULL
for (envir in names_envir_variables)
  {
    ok=FALSE
    while (!ok) {
      cat("Enter reaction norm model for variable",envir,"\n",
          "and demographic parameter", demographic_parameter,"\n",
          "(or 'h' for help) : ")
      shape <- append(shape,readline("Enter: ")) # prompt
      ok = (shape %in% allshapes)
      if (!ok) {cat("Please chose among: ","\n",
                   paste(allshapes,collapse="\n"))
                shape=shape[-length(shape)]
      }
    }
  }
shape
}

set_prior_vector_from_keyb <- function(name,n)
{
  # sets a vector of priors from keyboard
  # args: 
  # name of the vector to create
  # n: length of the vector
  #
ok = FALSE
while(!ok)
  {
  prior_dist <- readline(paste("Enter prior distribution for",
                 name,"(enter h for help): "))
  ok = prior_dist%in%c("uniform","log_uniform","normal","log_normal")
  if (!ok) cat("\n","models implemented are :",
               "\n","'uniform'",
               "\n","'log_uniform'",
               "\n","'normal'",
               "\n","'log_normal'"
               )  
  }
parameters_names_prior_dist=switch(prior_dist,
                           uniform = c("min","max"),
                           log_uniform = c("min","max"),
                           normal = c("mean","sd"),
                           log_normal = c("mean","sd"))
params = NULL
if (prior_dist=="log_normal") {cat("for log-normal, note that: ",
                                "\n"," mean and sd are on the log scale")}
if (prior_dist=="log_uniform") {cat("for log-uniform, note that: ",
                                 "\n","min and max are not on the log scale",
                                 "\n","but on the variable scale")}
for (paramname in parameters_names_prior_dist)
  {
  params = append(params,
                  readline(
                    paste("Enter ",paramname,
                        " for ",prior_dist," distribution: ")
                    )
                  )
  }
params=as.numeric(params)
names(params) = parameters_names_prior_dist
switch(prior_dist,
       uniform=runif(n,params["min"],params["max"]),
       log_uniform=exp(runif(n,log(params["min"]),log(params["max"]))),
       normal=rnorm(n,params["mean"],params["sd"]),
       log_normal=rlnorm(n,params["mean"],params["sd"])
         )
}

set_ref_table_from_keyb <- function(rasterStack,n)
{
  # function to create model parameter names of reference table
  # arg: none
  # value: reference table and average model 
  #
  # reaction norm models
  shapesK = input_reaction_norm_shape_model("K",names(rasterStack))
  pKnames = react_norm_param(shapesK)
  shapesr = input_reaction_norm_shape_model("r",names(rasterStack))
  prnames = react_norm_param(shapesr)
  # mutation model
  ok=FALSE
  while (!ok) {mutation_model <- readline("Enter mutation model (or 'h' for help) : ") # prompt
               ok = (mutation_model %in% c("tpm","bigeometric","stepwise"))
               if (!ok) {
                 cat("\n","models implemented are :",
                     "\n","'stepwise'",
                     "\n","'bigeometric'",
                     "\n","'tpm': two phase mutation model")
               }
  } 
  mut_param_names = switch(mutation_model,
                       bigeometric = "sigma2",
                       tmp = c("p","sigma2"),
                       stepwise = NULL
                       )
  # Dispersion
  ok=FALSE
  while (!ok) {shapeDisp <- readline("Enter dispersion model (or 'h' for help) : ") # prompt
               ok = (shapeDisp %in% c("fat_tail1","gaussian",
                                           "exponential","contiguous",
                                           "island", "fat_tail2",
                                           "gaussian_island_mix"))
               if (!ok) {
                 cat("\n","models implemented are :",
                     "\n","'fat_tail1' (Chapman et al)",
                     "\n","'gaussian'",
                     "\n","'exponential'",
                     "\n","'contiguous'",
                     "\n","'island'",
                     "\n","'fat_tail2' (Moilanen et al)",
                     "\n","'gaussian_island_mix'")
               }
  } 
  Dispersion_parameter_names = switch(shapeDisp,
                                      fat_tail1 =  c("alpha","beta"),
                                      gaussian = c("sd"),
                                      gaussian = c("sd"),
                                      exponential = c("mean"),
                                      contiguous = c("m"),
                                      island = c("m"),
                                      fat_tail2 = c("alpha","beta"),
                                      gaussian_island_mix = c("sd","m")
                                      )
# set priors
priors = list()
priors$shapesK = shapesK
priors$pK = rep(NA,length(pKnames));names(priors$pK)=pKnames
priors_names <- c(prnames,pKnames,mut_param_names,Dispersion_parameter_names)
df = as.data.frame(matrix(NA,nrow=n,ncol=length(priors_names)))
colnames(df)=priors_names
for (name in priors_names)
{
  df[,name] <- set_prior_vector_from_keyb(name,n)
}
df
}
  
input_priors <- function()
{
  # function to create prior values for reference table and average model
  # arg: none
  # value: reference table and average model 
  #
  nb_simul <- as.numeric(readline("Number of simulations: "))
  ok=FALSE
  while (!ok) {shape_model <- readline("Enter mutation model or 'h' for help : ") # prompt
               ok = (shape_model %in% c("tpm","bigeometric","stepwise"))
               if (!ok) {
                 cat("\n","models implemented are :",
                     "\n","'stepwise'",
                     "\n","'bigeometric'",
                     "\n","'tpm': two phase mutation model")
               }
  } 
  if (mutation_model=="bigeometric")
  {
    sigma2Dist <- readline("Enter distribution of variance maximum of geometric distribution (1/p): ")
    if (sigmaDist=="uniform")
    {
      sigma2Max <- readline("Enter variance maximum of geometric distribution (1/p): ")
      sigma2Min <- readline("Enter variance minimum of geometric distribution (1/p): ")
      sigma2 <- runif(n=nb_simul,min=sigma2Min,max=sigma2Max)
    }
  }
  if (mutation_model=="tpm")
  {
    sigma2Dist <- readline("Enter prior distribution shape for variance of geometric distribution (1/p): ")
    if (sigma2Dist=="uniform")
    {
      sigma2Max <- as.numeric(readline("Enter maximum of variance of geometric distribution: "))
      sigma2Min <- as.numeric(readline("Enter minimum of variance of geometric distribution: "))
      sigma2 <- runif(nb_simul,sigma2Min,sigma2Max)
      pMax <- readline("Enter stepwise maximum proportion: ")
      pMin <- readline("Enter stepwise minimum proportion: ")
      p <- runif(nb_simul,pMin,pMax)
    }
  }
}

set_model <- function(pK, pr, shapesK, shapesr, shapeDisp, pDisp,
                      mutation_rate, initial_genetic_value,
                      mutation_model,stepvalue,
                      mut_param)
{
  # sets a genetic and environemental demographic model as a list 
  # arg: parameters
  # value : list describing the models (shapes of distribution and parameters)
  model = list(pK=pK, pr=pr,
             shapesK=shapesK, shapesr=shapesr,
             shapeDisp=shapeDisp, pDisp=pDisp,
             mutation_rate=mutation_rate, 
             initial_genetic_value=initial_genetic_value, 
             mutation_model=mutation_model,stepvalue=stepvalue,
             mut_param=mut_param)
check_model(model)
}

check_dispersion_model <- function(shapeDisp,pDisp)
{
  # checks that shapeDisp and pDisp are compatible
  # arg: 
  # value:
  compatible = switch(shapeDisp,
                      fat_tail1 =  c("alpha","beta")%in%names(pDisp),
                      gaussian = c("sd")%in%names(pDisp),
                      gaussian = c("sd")%in%names(pDisp),
                      exponential = c("mean")%in%names(pDisp),
                      contiguous = c("m")%in%names(pDisp),
                      island = c("m")%in%names(pDisp),
                      fat_tail2 = c("alpha","beta")%in%names(pDisp),
                      gaussian_island_mix = c("sd","m")%in%names(pDisp)
  )
  if (!all(compatible)) {
    Message = switch(shapeDisp,
                     fat_tail1 =  "fat_tail1 : pDisp=c(alpha=..,beta=..)",
                     gaussian = "gaussian : pDisp=c(sd=..)",
                     exponential = "exponential: pDisp=c(mean=..)",
                     contiguous = "contiguous: pDisp=c(m=..)",
                     island = "island: pDisp=c(m=..)",
                     fat_tail2 = "fat_tail2: pDisp=c(alpha=..,beta=..)",
                     gaussian_island_mix = "gaussian_island_mix: pDisp=c(sd=..,m=..)"
    )
                     stop (paste("check dispersion model. For",Message))
  }
"Dispersion model OK"
}

check_reaction_norm <- function(shapes,p)
{
  # Checks whether shapes and parameter matrix of reaction norms
  # corresponds for all the variables
  #
  if (!all(shapes %in% c("enveloppe", "envelin", "envloglin","loG","linear",
                         "conquadratic","conquadraticskewed","conquadraticsq",
                         "conquadraticskewedsq","constant","linearPositive")))
    stop("shape of reaction norm is unknown, 
         please chose among 'constant', 'enveloppe', 'linearPositive' 'envelin', 'envloglin','loG','linear',
         'conquadratic','conquadraticskewed','conquadraticsq' or 'conquadraticskewedsq'
         ") else if (length(shapes)!=dim(p)[2])
           stop(paste("reaction norm shapes and parameters do not use the same number of
                       environemental variables.")) else {
                        compatible=NA
                        for (i in 1:length(shapes))
                        {
                          compatible[i] <- switch(shapes[i],
                                                  constant = all(c("Y")%in%row.names(p)),
                                                  enveloppe = all(c("Xmin","Xmax","Yopt")%in%row.names(p)),
                                                  envelin = all(c("Yxmin","Yxmax","Xmin","Xmax")%in%row.names(p)),
                                                  envloglin = all(c("Yxmin","Yxmax","Xmin","Xmax")%in%row.names(p)),
                                                  linear = all(c("X0","slope")%in%row.names(p)),
                                                  linearPositive = all(c("X0","slope")%in%row.names(p)),
                                                  conquadratic = all(c("Xmin","Xmax","Xopt")%in%row.names(p)),
                                                  conquadraticsq = all(c("Xmin","Xmax","Xopt")%in%row.names(p)),
                                                  conquadraticskewed = all(c("Xmin","Xmax","Xopt","Yopt")%in%row.names(p)),
                                                  conquadraticskewedsq = all(c("Xmin","Xmax","Xopt","Yopt")%in%row.names(p)))
                        } 
                        if (!all(compatible)) 
                        {
                        Message=NA
                        for (i in which(!compatible))
                          {
                          Message = switch(shapes[i],
                                           constant = "rownames of pDisp matrix parameter for constant shape is 'Y'. ",
                                           enveloppe = "rownames of pDisp matrix parameter for enveloppe shape are 'Xmin', 'Xmax' and 'Xopt'. ",
                                           envelin = "rownames of pDisp matrix parameter for envelin shape are 'Yxmin', 'Yxmax', 'Xmin', and 'Xmax'. ",
                                           envloglin = "rownames of pDisp matrix parameter for enveloglin shape are 'Yxmin', 'Yxmax', 'Xmin', and 'Xmax'. ",
                                           linear = "rownames of pDisp matrix parameter for linear shape are 'X0' and 'slope'. ",
                                           linearPositive = "rownames of pDisp matrix parameter for linearPositive shape are 'X0' and 'slope'. ",
                                           conquadratic = "rownames of pDisp matrix parameter for conquadratic shape are 'Xmin', 'Xmax' and Xopt. ",
                                           conquadraticsq = "rownames of pDisp matrix parameter for conquadraticsq shape are 'Xmin', 'Xmax' and Xopt. ",
                                           conquadraticskewed = "rownames of pDisp matrix parameter for conquadraticskewed shape are 'Xmin', 'Xmax', Xopt and Yxopt. ",
                                           conquadraticskewed = "rownames of pDisp matrix parameter for conquadraticskewed shape are 'Xmin', 'Xmax', Xopt and Yxopt. "
                          )
                          }
                        stop(paste("Error in reaction norm model settings:", 
                                   "for environmental variable(s) number",
                                   paste(which(!compatible),collapse=" and "),
                                   ", parameters corresponding to the shape annouced are not provided as rownames. Note that",
                                   paste(Message,collapse=" ")))
                        }
                      }
"reaction norm parameters OK"
}

check_mutation_model<-function(mutation_model,mutation_parameter)
{
  compatible = 
  switch(mutation_model,
         bigeometric = c("sigma2")%in%names(mutation_parameter),
         tmp = c("p","sigma2")%in%names(mutation_parameter)
         )
  if (!all(compatible)) {
    Message = switch(shapeDisp,
                     bigeometric =  "bigeometric : mut_param=c(sigma2=..)",
                     tmp = "tmp : mut_param=c(p=.., sigma2=..)"
                     )
    stop(paste("check mutation parameter(s) names for",Message))
  }
"mutation model OK"
}

check_model <- function(model)
{
check_reaction_norm(model$shapesr,model$pr)
check_reaction_norm(model$shapesK,model$pK)
check_mutation_model(model$mutation_model,model$mut_param)
check_dispersion_model(model$shapeDisp,model$pDisp)
}


setClass("DispersionModel", representation(ID="numeric",
                                           shapeDisp="character",
                                           pDisp="vector"
                                           )
         )


setClass("MutationModel",representation(model="character"
                                        )) 
setClass("EnvDemogenetModel",representation(ID = "numeric",
                                              K = "ReactionNorm",
                                              r = "ReactionNorm",
                                              Dispersion = "DispersionModel"
                                              ))





#
# add_genetic_to_coaltable function
# adds genetic values to a coalescent table containing mutation number per branch
# knowing initial genetic value of the ancastor and mutation model

genetics_of_coaltable <- function(coaltable,initial_genetic_value,mutation_model,stepvalue=2,mut_param=c(p=.5,sigma=2))
{
 switch(mutation_model,
        step_wise = stepwise(coaltable,initial_genetic_value,stepvalue),
        tpm = tpm(coaltable,initial_genetic_value,stepvalue,mut_param),
        bigeometric = bigeometric(coaltable,initial_genetic_value,stepvalue,mut_param)
        ) 
}

# coalescent_2_phylog
# converts a coalescent simulation to a class phylog tree (library ape)
#

coalescent_2_phylog <- function(coalescent)
{
read.tree(text=coalescent_2_newick(coalescent))
}

#
# plot_coalescent plots a coalescent simulation
# argument: output of simul_coalescent()

plot_coalescent <- function(coalescent,genetic_table,with_landscape=FALSE,rasK=NULL,legend_right_move=-.2)
{
  if (with_landscape) {par(mfrow=c(1,2),oma=c(0,0,0,4),xpd=TRUE)}else{par(mfrow=c(1,1),oma=c(0,0,0,4),xpd=TRUE)}
  tipcells <- geneticData$Cell_numbers[as.numeric(coalescent_2_phylog(coalescent)$tip.label)]
  tipcols = rainbow(ncell(rasK))[tipcells]
  phylog_format_tree <- coalescent_2_phylog(coalescent)
  phylog_format_tree$tip.label <- paste(phylog_format_tree$tip.label,genetic_table[order(genetic_table$coalescing)[as.numeric(phylog_format_tree$tip.label)],"genetic_value"],sep=":")
  plot(phylog_format_tree,direction="downward",tip.color=tipcols)
  legend("topright", title="demes", cex=0.75, pch=16, col=tipcols[!duplicated(tipcols)], legend=tipcells[!duplicated(tipcols)], ncol=2, inset=c(legend_right_move,0))
  if (with_landscape) {plot(rasK)}
}

# summary_stat calculates summary stats for observed and simulated data and creates a reference table
# geneticDataSimulList : a list of simulations, with sublist geneticData and sublist log_lik_forward 
#


set_priors <- function(variables,Min,Max,nb_lines)
{
  df = data.frame(matrix(NA,nrow=1,ncol=length(variables)))
  for (i in variables)
  {
    
  }
}

new_reference_table <- function(geneticData,Distance,priors)
{
  # creates a new reference table from genetic data and priors data
  # !!!!!! details : first line of reference table is NOT rotated observed genetic data
  #list of of parameters
  
  #list of summary stats
  GenetDist = switch(Distance,
                     Goldstein = dist(geneticData[,grep("Locus",colnames(geneticData))])^.5,
                     pID = pID(geneticData)
  )
  Significant_Components <- which(abs(prcomp(GenetDist)$x)>1E-5)
  df <- as.data.frame(matrix(NA,nrow=1,ncol=length(Significant_Components)+1,dimnames=list("l1",c(paste("C",c(Significant_Components),sep=""),"loglik"))))[FALSE,]
  df
}


add_summary_stat <- function(reference_table,geneticDataSim,rotation,forward_log_lik,Distance="Goldstein")
{

  # add_summary_stats
  # add summary statistics of a simulation to reference table for ABC analysis
  # arguments: 
  # reference_table the reference table to append
  # rotation tha PCA rotation to apply to geneticData
  # geneticDataSim : the simulated genetic data
  # forward_log_lik : forward log likelihood of the simulated genetic data
  # Distance : distance to apply to genetic data (Goldstein : difference in number of repeats)
  # pID : proportion of identity.
  # 
  # value: appended reference_table
  
  #1) We calculate a matrix of observed genetic distance between individuals (DSAI = sum_j (Shared alleles) / (number of loci))
  #   or reapeat number distance (Goldstein 1995)
  GenetDist = switch(Distance,
                     Goldstein = dist(geneticData[,grep("Locus",colnames(geneticData))])^.5,
                     pID = pID(geneticData)
  )
  #2) We express simulated data in these axis. First summary stats are the values of each invidiual
  # in each of these major axes
  summary_stats = as.matrix(GenetDist) %*% rotation
rbind(reference_table,as.vector(summary_stats))
}

fill_reference_table <- function(geneticData,Distance,rasterStack=rasterStack,
                                 pK=pK, pr=pr,
                                 shapesK=shapesK, shapesr=shapesr,
                                 shapeDisp=shapeDisp, pDisp=pDisp,
                                 mutation_rate=1E-1, 
                                 initial_genetic_value=initial_genetic_value, 
                                 mutation_model="tpm",stepvalue=2,
                                 mut_param=c(p=.5,sigma2=4))
{
  # filling a reference table using parameters values
  # args:
  # parameters of the model, geneticData, type of distance used
  rotation = PCA_rotation(geneticData)
  ref_table = new_reference_table(geneticData,Distance="Goldstein")
  simulated_genetic <- geneticData[,grep("locus",colnames(geneticData))]
  # simulation de genetic data pour chaque locus
  for (i in 1:length(grep("Locus",colnames(geneticData))))
  {
    new_simulation <- simul_coalescent(geneticData=geneticData,
                                       rasterStack=rasterStack,
                                       pK=pK, pr=pr,
                                       shapesK=shapesK, shapesr=shapesr,
                                       shapeDisp=shapeDisp, pDisp=pDisp,
                                       mutation_rate=1E-1, 
                                       initial_genetic_value=initial_genetic_value, 
                                       mutation_model="tpm",stepvalue=2,
                                       mut_param=c(p=.5,sigma2=4))
    simulated_genetic[,i] <- new_simulation$genetic_values[order( new_simulation$coalescing)[1:dim(geneticData)[1]],"genetic_value"]
    forward_log_lik <- new_simulation$forward_log_prob
  }
  simulated_genetic  
  ref_table = add_summary_stat(reference_table=ref_table,
                               geneticDataSim=simulated_genetic,
                               rotation=rotation,
                               forward_log_lik=,Distance="Goldstein")
  ref_table
  
}



#
# genetic_simulation
# this function simulates genetic data from a coalescent
#
#


test_stabilite_a_value <- function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice){
  ## ref: Rousset et al. J Evol Biol,13  (2000) 58-62.
  vecteur_a_value <-c(0)
  GenDist=list()
  for(i in 1: nb_generations){#i=1
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
}

# to get a matrix of a-values between all cells of a landscape from 
# simulations of evolution where not all cells necesarily contain individuals
#
#

a_value_matrix_from_forward_simulation <- function(geneticFinal,rasterStack)
{
  matrixQb = (1-Qwithin_pair(geneticFinal)+2*(Qwithin_pair(geneticFinal)-Qbetween(geneticFinal,dim(geneticFinal))))
  matrixQw = 2*(1-Qwithin_pop(geneticFinal))
  vecteur_a_value = matrixQb/matrixQw-1/2
  vecteur_a_value[is.na(vecteur_a_value)] <-0
  cell_numbers_list <- levels(as.factor(geneticFinal$Cell_numbers))
  agg <- aggregate(vecteur_a_value,by=list(geneticFinal$Cell_numbers),FUN="mean")[,-1]
  agg2 <- t(aggregate(t(agg),by=list(geneticFinal$Cell_numbers),FUN="mean"))[-1,]
  row.names(agg2) <- as.numeric(cell_numbers_list)
  colnames(agg2) <- as.numeric(cell_numbers_list)
  agg3 <- matrix(NA, nrow=ncell(rasterStack),ncol=ncell(rasterStack))
  agg3[as.numeric(row.names(agg2)),as.numeric(colnames(agg2))] <- agg2
agg3
}

# Function that determine number of dispersal parameters from dispersal shapeDisp used
# Useful for nlm estimation
nbpDisp <- function(shapeDisp){
  (switch(shapeDisp,
          fat_tail1 = 2,
          gaussian = 1,
          exponential = 1,
          contiguous = 1,
          island = 1,
          fat_tail2 = 2))
}

# Function that simulate a genetic data with parameters given. 
#It returns a list of 2 variables:  final genetic data observed and matrix of a-value observed
simulationGenet <- function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, mutationRate, nbLocus, initial_locus_value, shapeDisp, pDisp, nb_generations,ind_per_cell=30){
  K <- subset(ReactNorm(values(donneesEnvironmentObs), pK , shapesK),select="Y") # carrying capacity
  r <- subset(ReactNorm(values(donneesEnvironmentObs), pR , shapesR),select="Y") # growth rate
  Rast_K <- donneesEnvironmentObs ; values(Rast_K) <- K
  Rast_r <- donneesEnvironmentObs ; values(Rast_r) <- r
  geneticData = CreateGenetArray(Rast_K, nbLocus,initial_locus_value,Option="full_pop")
  geneticData[,"Cell_number_init"] <- geneticData[,"Cell_numbers"]
  dimGeneticData = dim(geneticData)
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
  transitionmatrice = transitionMatrixBackward(Npop = K, migration= migrationM)
  geneticDataFinal = test_stabilite_a_value(geneticData, mutationRate, dimGeneticData, nb_generations,transitionmatrice)
  a_value_obs = geneticDataFinal[[2]]
  geneticDataFinal = geneticDataFinal[[1]]
  return(list(geneticDataFinal, a_value_obs))
}

# plots matrixes of forward mutation for validation
#
#

matrixes_forward <- function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, shapeDisp, pDisp, a_value_obs, a_value_att, file=NULL, mutationRate,nbLocus, initial_locus_value,indpercell)
{
  nblayers =dim(donneesEnvironmentObs)[3]
  nCell = ncell(donneesEnvironmentObs)
  Cell_numbers <- 1:nCell
  geneticObs = simulationGenet(donneesEnvironmentObs,pK,pR,shapesK,shapesR,mutationRate,nbLocus,initial_locus_value,shapeDisp,pDisp,nb_generations=5000,indpercel)
  finalGenetData = geneticObs[[1]]
  a_value_simul = geneticObs[[2]]
  K <-reactNorm(donneesEnvironmentObs, pK, shapesK) # carrying capacity
  r <-reactNorm(donneesEnvironmentObs, pR, shapesR) # carrying capacity
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
  transitionmatrice = transitionMatrixForward(r, K, migration= migrationM)
  land_size <- raster(matrix(K[1,],nrow=dim(donneesEnvironmentObs)[1],ncol=dim(donneesEnvironmentObs)[2]))# be careful transposed
  land_size <- raster(matrix(K[1,]))
  extent(land_size) <- extent(donneesEnvironmentObs)
  list(land_size,transitionmatrice,a_value_simul,a_value_att)
}

# MAIN DES FONCTIONS POUR SIMULER DES DONNEES GENETIQUES OBSERVEES
#mutationRate = 10^(-4);                             shapeDisp = "fat_tail1";
#nbLocus=10;                                         initial_locus_value=200
#pDisp = c(alphaDisp = 20, betaDisp = 1.5);     nbpDisp =nbpDisp(shapeDisp);
#aggregate_factor=16
#donneesEnvironmentObs = Aggregate_and_adjust_raster_to_data(raster(paste(wd,envdir,envfiles,sep="")),release=read.table(paste(wd,genetfile,sep="")), recovery=read.table(paste(wd,genetfile,sep="")), extend_band_size=1, aggregate_index=aggregate_factor)
  
# donneesEnvironmentObs = raster(matrix(NA,nrow=20,ncol=20))
#values(donneesEnvironmentObs) <-log(ceiling(runif(ncell(donneesEnvironmentObs),0,1)*3))
#plot(donneesEnvironmentObs)
#nblayers =dim(donneesEnvironmentObs)[3]
#nCell = ncell(donneesEnvironmentObs)
#geneticObs = simulationGenet(donneesEnvironmentObs,alpha=0, beta =1,mutationRate,nbLocus, initial_locus_value,shapeDisp,pDisp,nb_generations=5000)
#finalGenetData = geneticObs[[1]]
#a_value_obs = geneticObs[[2]]

#####################################
##### ESTIMATION DES PARAMETRES #####
#####################################
expect_a_value <- function(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,Cell_numbers,nbpDisp,nblayers,shapeDisp)
{
  K = ReactNorm(rasterStack =donneesEnvironmentObs, pK, shapesK)
  r = ReactNorm(rasterStack =donneesEnvironmentObs, pR, shapesR)
  matrice_migration = migrationMatrix(donneesEnvironmentObs,shapeDisp, shapeDisp)
  matrice_transition = transitionMatrixBackward(K,r, matrice_migration)
  matrice_laplace = laplaceMatrix(matrice_transition)
  commute_time = resistDist(matrice_laplace)
  genetic_dist_att = geneticDist(commute_time, popSize)
  a_value_att = genetic_dist_att[Cell_numbers,Cell_numbers]
  a_value_att
}

ssr <- function(p){
  a_value_att = expect_a_value(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,Cell_numbers=finalGenetData$Cell_Number,nbpDisp,nblayers,shapeDisp)
  return(mean((a_value_obs - a_value_att)^2))
}

#initial = c(0,1,20,1.5)
#fct_erreur_calc = nlm(f = ssr , p = initial, hessian = FALSE)
#p=initial
#a_value_graph_model = expect_a_value(donneesEnvironmentObs,initial,finalGenetData$Cell_Number,nbpDisp=2,nblayers=1,shapeDisp="fat_tail1")
#ssr(initial)
#parametres_infere = fct_erreur_calc$estimate
#parametres_reels = c(alpha = 0, beta = 2,beta = 1, alphaDisp = 20, betaDisp = 1.5)

validation <- function(donneesEnvironmentObs,pK = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))), pR = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))),shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"),shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"), shapeDisp="fat_tail1", pDisp = c(mean=0.32,shape=1.6), file=NULL,mutationRate,nbLocus, initial_locus_value,nb_generations=5000,indpercell=30)
{
  nblayers =dim(donneesEnvironmentObs)[3]
  nCell = ncell(donneesEnvironmentObs)
  Cell_numbers <- 1:nCell
  K <- ReactNorm(values(donneesEnvironmentObs), pK , shapesK)$Y # recuperation de la carrying capacity
  r <- ReactNorm(values(donneesEnvironmentObs), pR , shapesR)$Y # recuperation du growth rate
  # Migration function used to calculate descendant from parents (draw with replacement)
  migrationM = migrationMatrix(donneesEnvironmentObs,shapeDisp,pDisp)
  transitionmatrice = transitionMatrixBackward(Npop = K, migration= migrationM)
  geneticObs = simulationGenet(donneesEnvironmentObs,alpha, beta, mutationRate,nbLocus, initial_locus_value,shapeDisp,pDisp,nb_generations=5000,indpercell)
  finalGenetData = geneticObs[[1]]
  land_size <- raster(matrix(K[1,],nrow=dim(donneesEnvironmentObs)[1],ncol=dim(donneesEnvironmentObs)[2]))
  a_value_simul = a_value_matrix_from_forward_simulation(geneticObs[[1]],land_size)
  a_value_theory_stepping_stone_model = distanceMatrix(donneesEnvironmentObs)/(4*0.05)
  a_value_theory_island_model <- matrix(pDisp[1],nrow=nCell,ncol=nCell)-pDisp[1]*diag(nCell)
  a_value_theory_graph_model <- expect_a_value(donneesEnvironmentObs,pK,pR,shapesK,shapesR,pDisp,Cell_numbers,nbpDisp=nbpDisp,nblayers=nblayers,shapeDisp)
  if (!is.null(file)){jpeg(file)}
  par(mfrow=c(2,3))
  plot(land_size,main="Population sizes")
  #plot(raster(migrationM))
  plot(raster(transitionmatrice),main="Transition matrix")
  plot(raster(a_value_simul),main="Simul genet differ")
  plot(raster(a_value_theory_stepping_stone_model),main="Expected stepping stone")
  plot(raster(a_value_theory_island_model),main="Expected island")
  plot(raster(a_value_theory_graph_model),main="Expected graph model")
  if (!is.null(file)){dev.off()}
  list(land_size,transitionmatrice,a_value_simul,a_value_theory_stepping_stone_model,a_value_theory_island_model,a_value_theory_graph_model,finalGenetData)
}


#test nlm
ssr <- function(p){
  popSize = K_Function(rasterStack = donneesEnvironmentObs, p[1], p[2:(nblayers+1)])
  matrice_migration = migrationMatrix(donneesEnvironmentObs,shapeDisp, p[(nblayers+2):(nbpDisp+nblayers+1)])
  matrice_transition = transitionMatrixBackward(popSize, matrice_migration)
  matrice_laplace = laplaceMatrix(matrice_transition)
  commute_time = resistDist(matrice_laplace)
  genetic_dist_att = geneticDist(commute_time, p[(nbpDisp+nblayers+2)])
  
  a_value_att = genetic_dist_att[finalGenetData$Cell_Number,finalGenetData$Cell_Number]
  return(mean((a_value_obs - a_value_att)^2))
}



