setGeneric(
  name = "laplaceMatrix",
  def=function(object){return(standardGeneric("laplaceMatrix"))}
)


setGeneric(
  name = "ordinary_laplacian",
  def=function(object){return(standardGeneric("ordinary_laplacian"))}
)

setGeneric(
  name = "commute_time_undigraph",
  def=function(object){return(standardGeneric("commute_time_undigraph"))}
)

setGeneric(
  name = "hitting_time_digraph",
  def=function(object){return(standardGeneric("hitting_time_digraph"))}
)

setGeneric(
  name = "genetDistUndigraph",
  def=function(object,popSize,mutation_rate){return(standardGeneric("genetDistUndigraph"))}
)

setGeneric(
  name = "genetDistDigraph",
  def=function(object,popSize,mutation_rate,method="Goldstein95"){return(standardGeneric("genetDistDigraph"))}
)

#############  A MODIFIER   #####################


setGeneric(
  name = "a_matrix",
  def=function(gDat){return(standardGeneric("a_matrix"))}
)


setGeneric(
  name = "absorbingTransitionComplete",
  def=function(transition,N){return(standardGeneric("absorbingTransitionComplete"))}
)


setGeneric(
  name = "coalescenceTimeProbDistrib",
  def=function(Qlist){return(standardGeneric("coalescenceTimeProbDistrib"))}
)


setGeneric(
  name = "enveloppe",
  def=function(X,p){return(standardGeneric("enveloppe"))}
)


setGeneric(
  name = "add_br_length_and_mutation",
  def=function(coalescent,mutation_rate,initial_locus_value,allelenames){return(standardGeneric("add_br_length_and_mutation"))}
)


setGeneric(
  name = "absorbingTransition",
  def=function(transition,N){return(standardGeneric("absorbingTransition"))}
)


setGeneric(
  name = "summary_stat",
  def=function(geneticDataObs,geneticDataSimulList,log_lik_simul_list){return(standardGeneric("summary_stat"))}
)


setGeneric(
  name = "coalescent_2_newick",
  def=function(coalescent){return(standardGeneric("coalescent_2_newick"))}
)


setGeneric(
  name = "K_Function",
  def=function(rasterStack, p, shapes){return(standardGeneric("K_Function"))}
)


setGeneric(
  name = "ordinary_laplacian",
  def=function(transition){return(standardGeneric("ordinary_laplacian"))}
)


setGeneric(
  name = "conquadratic",
  def=function(X,p){return(standardGeneric("conquadratic"))}
)


setGeneric(
  name = "Qbetween",
  def=function(gDat){return(standardGeneric("Qbetween"))}
)


setGeneric(
  name = "probgenet",
  def=function(transition,gDat){return(standardGeneric("probgenet"))}
)


setGeneric(
  name = "gridRepnDispFunction",
  def=function(dynamics,r,K,d=.9,ptG, migration,overlapping=TRUE){return(standardGeneric("gridRepnDispFunction"))}
)


setGeneric(
  name = "Genetic_Dist",
  def=function(gDat,method="Goldstein"){return(standardGeneric("Genetic_Dist"))}
)


setGeneric(
  name = "migrationMatrix",
  def=function(rasterStack,shapeDisp, pDisp){return(standardGeneric("migrationMatrix"))}
)


setGeneric(
  name = "R_Function",
  def=function(rasterStack, alpha, beta){return(standardGeneric("R_Function"))}
)


setGeneric(
  name = "genetics_of_coaltable",
  def=function(coaltable,initial_locus_value,mutation_model="stepwise",stepvalue=2){return(standardGeneric("genetics_of_coaltable"))}
)


setGeneric(
  name = "OneCol23Dims",
  def=function(gDat){return(standardGeneric("OneCol23Dims"))}
)


setGeneric(
  name = "fundamentalMatrixAbsorbingChain",
  def=function(transientTransitionMatrix){return(standardGeneric("fundamentalMatrixAbsorbingChain"))}
)


setGeneric(
  name = "commute_time_undigraph",
  def=function(matrice_transition){return(standardGeneric("commute_time_undigraph"))}
)


setGeneric(
  name = "linearizedFstDigraph1",
  def=function(transition, popSize){return(standardGeneric("linearizedFstDigraph1"))}
)


setGeneric(
  name = "coalist_2_coaltable",
  def=function(coalist,locusnames,initial_locus_value){return(standardGeneric("coalist_2_coaltable"))}
)


setGeneric(
  name = "comuteTimeDigraph",
  def=function(transition, popSize){return(standardGeneric("comuteTimeDigraph"))}
)


setGeneric(
  name = "a_value_matrix",
  def=function(gDat,rasterStack){return(standardGeneric("a_value_matrix"))}
)


setGeneric(
  name = "genetDistAbsorbingMethod",
  def=function(transition,N,mutation_rate){return(standardGeneric("genetDistAbsorbingMethod"))}
)


setGeneric(
  name = "plot_coalescent",
  def=function(coalescent_simulated,with_landscape=FALSE,rasK=NULL,legend_right_move=-.2,file=NA){return(standardGeneric("plot_coalescent"))}
)


setGeneric(
  name = "nbpDisp",
  def=function(shapeDisp){return(standardGeneric("nbpDisp"))}
)


setGeneric(
  name = "envelin0",
  def=function(X,p,log=FALSE){return(standardGeneric("envelin0"))}
)


setGeneric(
  name = "TwoCols2OneCol",
  def=function(tip_genotype){return(standardGeneric("TwoCols2OneCol"))}
)


setGeneric(
  name = "Qbetween_2",
  def=function(gDat){return(standardGeneric("Qbetween_2"))}
)


setGeneric(
  name = "SSw",
  def=function(gDat){return(standardGeneric("SSw"))}
)


setGeneric(
  name = "distanceMatrix",
  def=function(rasterStack){return(standardGeneric("distanceMatrix"))}
)


setGeneric(
  name = "transitionMatrixForward",
  def=function(r,K, migration, meth="non_overlap"){return(standardGeneric("transitionMatrixForward"))}
)


setGeneric(
  name = "genetDistDigraph",
  def=function(transition,popSize,mutation_rate,method="Goldstein95"){return(standardGeneric("genetDistDigraph"))}
)


setGeneric(
  name = "stepwise",
  def=function(coaltable,initial_locus_value,stepvalue,locusnames){return(standardGeneric("stepwise"))}
)


setGeneric(
  name = "coalescenceProb",
  def=function(Qlist,time){return(standardGeneric("coalescenceProb"))}
)


setGeneric(
  name = "Aggregate_and_adjust_raster_to_data",
  def=function(Envir_raster_stack,release,recovery,extend_band_size,aggregate_index){return(standardGeneric("Aggregate_and_adjust_raster_to_data"))}
)


setGeneric(
  name = "plot_genealogy",
  def=function(genealogy,file=NA){return(standardGeneric("plot_genealogy"))}
)


setGeneric(
  name = "OneCol2TwoCols",
  def=function(gDat){return(standardGeneric("OneCol2TwoCols"))}
)


setGeneric(
  name = "get_nj_tree",
  def=function(tip_genotype,mutation_model,step_value){return(standardGeneric("get_nj_tree"))}
)


setGeneric(
  name = "mmute",
  def=function(cells=c(1,2),transitionmatrice){return(standardGeneric("mmute"))}
)


setGeneric(
  name = "repnDispMutFunction",
  def=function(geneticData, dimGeneticData, mutationRate, transitionmatrice){return(standardGeneric("repnDispMutFunction"))}
)


setGeneric(
  name = "simul_coalescent",
  def=function(transitionList,geneticData){return(standardGeneric("simul_coalescent"))}
)


setGeneric(
  name = "genetDistUndigraph",
  def=function(transition,popSize,mutation_rate){return(standardGeneric("genetDistUndigraph"))}
)


setGeneric(
  name = "ocur",
  def=function(cells=c(1,2),transitionmatrice){return(standardGeneric("ocur"))}
)


setGeneric(
  name = "genetDist",
  def=function(tip_genotype,method="deltaMu",stepvalue=2,byCell=TRUE){return(standardGeneric("genetDist"))}
)


setGeneric(
  name = "ReactNorm",
  def=function(X,p,shapes){return(standardGeneric("ReactNorm"))}
)


setGeneric(
  name = "matrixes_forward",
  def=function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, shapeDisp, pDisp, a_value_obs, a_value_att, file=NULL, mutationRate,nbLocus, initial_locus_value,indpercell){return(standardGeneric("matrixes_forward"))}
)


setGeneric(
  name = "checkTimeInterval",
  def=function(transition = result$transitionmatrice,rasK = rasK,threshold = 1E-5){return(standardGeneric("checkTimeInterval"))}
)


setGeneric(
  name = "Qwithin_pop",
  def=function(gDat){return(standardGeneric("Qwithin_pop"))}
)


setGeneric(
  name = "proportional",
  def=function(X,p,Log=FALSE){return(standardGeneric("proportional"))}
)


setGeneric(
  name = "samplePrior",
  def=function(prior,method="random"){return(standardGeneric("samplePrior"))}
)


setGeneric(
  name = "formatGeneticData",
  def=function(gDat,rasK){return(standardGeneric("formatGeneticData"))}
)


setGeneric(
  name = "laplaceMatrix",
  def=function(transitionMatrix){return(standardGeneric("laplaceMatrix"))}
)


setGeneric(
  name = "validation",
  def=function(donneesEnvironmentObs,pK = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))), pR = matrix(c(100,400,300,10,0.5,1,300,2000,1500,10,.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12"))),shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"),shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed"), shapeDisp="fat_tail1", pDisp = c(mean=0.32,shape=1.6), file=NULL,mutationRate,nbLocus, initial_locus_value,nb_generations=5000,indpercell=30){return(standardGeneric("validation"))}
)


setGeneric(
  name = "Qwithin_pair",
  def=function(gDat){return(standardGeneric("Qwithin_pair"))}
)


setGeneric(
  name = "expect_linearizedFst_undigraph",
  def=function(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA){return(standardGeneric("expect_linearizedFst_undigraph"))}
)


setGeneric(
  name = "transitionMatrice",
  def=function(rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp){return(standardGeneric("transitionMatrice"))}
)


setGeneric(
  name = "degree2km",
  def=function(rasterStack){return(standardGeneric("degree2km"))}
)


setGeneric(
  name = "conquadraticskewedsq",
  def=function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12")))){return(standardGeneric("conquadraticskewedsq"))}
)


setGeneric(
  name = "populationSize",
  def=function(donneesEnvironmentObs, p, shapes){return(standardGeneric("populationSize"))}
)


setGeneric(
  name = "a_value_ind",
  def=function(gDat){return(standardGeneric("a_value_ind"))}
)


setGeneric(
  name = "envelinear",
  def=function(X,p,log=FALSE){return(standardGeneric("envelinear"))}
)


setGeneric(
  name = "ssr",
  def=function(p){return(standardGeneric("ssr"))}
)


setGeneric(
  name = "simulationGenet",
  def=function(donneesEnvironmentObs, pK, pR, shapesK, shapesR, mutationRate, nbLocus, initial_locus_value, shapeDisp, pDisp, nb_generations,ind_per_cell=30){return(standardGeneric("simulationGenet"))}
)


setGeneric(
  name = "coalescence_prob_time_distribution_matrix",
  def=function(transition,max_time_interval=4,rasK,threshold=1E-6){return(standardGeneric("coalescence_prob_time_distribution_matrix"))}
)


setGeneric(
  name = "linearizedFstUndigraph",
  def=function(transition, popSize){return(standardGeneric("linearizedFstUndigraph"))}
)


setGeneric(
  name = "plot_coal_time_depending_on_time_interval",
  def=function(meanCoalTimes,filen=NA,corr=TRUE,land=TRUE){return(standardGeneric("plot_coal_time_depending_on_time_interval"))}
)


setGeneric(
  name = "conquadraticsq",
  def=function(X,p){return(standardGeneric("conquadraticsq"))}
)


setGeneric(
  name = "Qwithin_pair_2",
  def=function(gDat){return(standardGeneric("Qwithin_pair_2"))}
)


setGeneric(
  name = "deltaMu",
  def=function(tip_genotype,stepvalue=2){return(standardGeneric("deltaMu"))}
)


setGeneric(
  name = "TwoCols2OneCol_",
  def=function(gDat){return(standardGeneric("TwoCols2OneCol_"))}
)


setGeneric(
  name = "conquadraticskewed",
  def=function(X,p=matrix(c(100,400,250,1,0.5,1,300,2000,1000,1,0.5,1),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yopt","Yxmin","Yxmax"),c("BIO1","BIO12")))){return(standardGeneric("conquadraticskewed"))}
)


setGeneric(
  name = "genetDistUndigraphForNLM",
  def=function(parameters){return(standardGeneric("genetDistUndigraphForNLM"))}
)


setGeneric(
  name = "expect_linearizedFst_digraph",
  def=function(rasterStack,pK,pR,shapesK,shapesR,pDisp,shapeDisp,Cell_numbers=NA){return(standardGeneric("expect_linearizedFst_digraph"))}
)


setGeneric(
  name = "FstatsRas",
  def=function(gDat,by="cell",all=TRUE,cells=NULL){return(standardGeneric("FstatsRas"))}
)


setGeneric(
  name = "test_stabilite_a_value",
  def=function(geneticData, mutationRate, dimGeneticData, nb_generations=5000,transitionmatrice){return(standardGeneric("test_stabilite_a_value"))}
)


setGeneric(
  name = "likelihoodCoalescent",
  def=function(coalescent){return(standardGeneric("likelihoodCoalescent"))}
)


setGeneric(
  name = "SharedAlleleDistance",
  def=function(tip_genotype){return(standardGeneric("SharedAlleleDistance"))}
)


setGeneric(
  name = "SSb",
  def=function(gDat){return(standardGeneric("SSb"))}
)


setGeneric(
  name = "absorbingTime",
  def=function(fundamentalMatrixAbsChain,Qline){return(standardGeneric("absorbingTime"))}
)


setGeneric(
  name = "S",
  def=function(rasterStack){return(standardGeneric("S"))}
)


setGeneric(
  name = "id_mat_prob",
  def=function(A,B){return(standardGeneric("id_mat_prob"))}
)


setGeneric(
  name = "tmra",
  def=function(coalescent_simulated){return(standardGeneric("tmra"))}
)


setGeneric(
  name = "CreateGenetArray",
  def=function(rasK, nb_locus, initial_locus_value,Option="sample_1_col_diploid",nind="Ne"){return(standardGeneric("CreateGenetArray"))}
)


setGeneric(
  name = "coalescent_2_phylog",
  def=function(coalescent){return(standardGeneric("coalescent_2_phylog"))}
)


setGeneric(
  name = "deme_coocurence_probability",
  def=function(pGenes,transition,time){return(standardGeneric("deme_coocurence_probability"))}
)


setGeneric(
  name = "linearizedFstDigraph",
  def=function(transition, popSize){return(standardGeneric("linearizedFstDigraph"))}
)


setGeneric(
  name = "linear",
  def=function(X,p,Log=FALSE){return(standardGeneric("linear"))}
)


setGeneric(
  name = "hitting_time_digraph",
  def=function(transition){return(standardGeneric("hitting_time_digraph"))}
)


setGeneric(
  name = "get_ultrametric_nj_tree",
  def=function(tip_genotype,mutation_model,step_value){return(standardGeneric("get_ultrametric_nj_tree"))}
)


setGeneric(
  name = "genetDistStepping",
  def=function(migration,popSize,mutation_rate){return(standardGeneric("genetDistStepping"))}
)


setGeneric(
  name = "forward_simul_landpopsize",
  def=function(N0,p, migration){return(standardGeneric("forward_simul_landpopsize"))}
)


setGeneric(
  name = "landGenetAnalysis",
  def=function(genetData,environmentalData,priors){return(standardGeneric("landGenetAnalysis"))}
  )


setGeneric(
  name = "transitionMatrixBackward",
  def=function(r,K, migration){return(standardGeneric("transitionMatrixBackward"))}
)


setGeneric(
  name = "simul_coal_genet",
  def=function(geneticData,rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp,mutation_rate=1E-2,initial_locus_value=200,mutation_model="stepwise",stepvalue=2,locusnames=NA){return(standardGeneric("simul_coal_genet"))}
)

