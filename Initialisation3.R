rm(list=ls())
wd1="/home/dupas/Bureau/" # fixe
wd2="/home/dupas/Bureau/Graph_Pop" # portable
wd3="/media/2To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # portable
wd4="/media/dupas/2To/IRD/ARTICLES/Dupas/NicheCoal/Graph_Pop" # fixe
source.file=c("graphPOP_0.117.R")
if (source.file %in% list.files(wd1)) {wd=wd1 ;setwd(wd);source("classes.R");source("generics.R");source("methods.R");source(source.file)}
if (source.file %in% list.files(wd2)) {wd=wd2 ;setwd(wd);source("classes.R");source("generics.R");source("methods.R");source(source.file)}
if (source.file %in% list.files(wd3)) {wd=wd3 ;setwd(wd);source("classes.R");source("generics.R");source("methods.R");source(source.file)}
if (source.file %in% list.files(wd4)) {wd=wd4 ;setwd(wd);source("classes.R");source("generics.R");source("methods.R");source(source.file)}
load("result20150127.Rdata")

# test simulation backward
# for each demographic varaibles (r and K)
# we have a niche function (reaction norm)
# the shapes of the reaction =norms are given by the variables shapeK and shaper
# the parameters of the reaction norm are given by  and pR

# we also have a shape of the dispersion function given here by the variable "shapeDisp"
# and dispersion parameters of the dispersion function
#pDisp = c(1/19,1)
#shapeDisp="fat_tail1"
mutation_model="stepwise";locusnames=NA;stepvalue=2;shapeDisp="contiguous";pDisp = .05;
mutation_rate=1E-2;nind="Ne";Option="sample_1col_diploid"
filen=NULL;mutationRate=1E-2;nbLocus=6; initial_locus_value=200

N=150;r=10
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
#shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
#pR = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
shapesR=c(BIO1="constant",BIO12="constant")
pR = matrix(c(r,r),nrow=1,ncol=2,dimnames=list(c("Y"),c("BIO1","BIO12")))



#Data2 <- data.frame(BIO1=sample(300:101,15),BIO12=sample(31:230*10,15));dim(Data2)
Data2 <- data.frame(BIO1=rep(120,10),BIO12=900,10);dim(Data2)
steppingStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=dim(Data2)[1]),xmn=0,xmx=dim(Data2)[1],ymn=0,ymx=1),"BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=dim(Data2)[1]),xmn=0,xmx=dim(Data2)[1],ymn=0,ymx=1)))
rasK=steppingStack;values(rasK)= as.matrix(ReactNorm(values(steppingStack),pK,shapesK)[,"Y"])
plot(rasK)
rasK
class(rasK)

# rasterStack=rasterStack;pK = pK; pR = pR;shapesK=shapesK;shapesR=shapesR; shapeDisp="contiguous"; pDisp = .5;filen="sinusline20.jpg";mutationRate=2E-2;nbLocus=12; stepvalue=2;initial_locus_value=200;mutation_model="stepwise";mutation_rate=2E-2;nind="Ne";locusnames=NA
result=validation_with_coalescent(rasterStack=steppingStack,pK = pK, pR = pR,shapesK=shapesK,
                                  shapesR=shapesR, shapeDisp="contiguous", pDisp = .1, 
                                  filen="stepping10.jpg",nbLocus=24, 
                                  stepvalue=2,initial_locus_value=200,mutation_model="stepwise",
                                  mutation_rate=2E-2,nind="Ne",stat="F")

Data2 <- data.frame(BIO1=310+90*sin(pi*(0:9)/2.5),BIO12=1400+950*sin(pi*(0:9)/2.5));dim(Data2)
N=50
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
sinusStepStack <- stack(raster(matrix(310+90*sin(pi*(0:32)/10),ncol=33),xmn=0,xmx=33,ymn=0,ymx=1),
                        raster(matrix(1400+950*sin(pi*(0:32)/10),ncol=33),xmn=0,xmx=33,ymn=0,ymx=1))
names(sinusStepStack)<-c("BIO1","BIO12")
rasK=sinusStepStack;values(rasK)= as.matrix(ReactNorm(values(sinusStepStack),pK,shapesK)[,"Y"])
plot(rasK)
# rasterStack=sinusStepStack;pK = pK; pR = pR;shapesK=shapesK; shapesR=shapesR; shapeDisp="contiguous_long_dist_mixt"; pDisp = c(plongdist=.05,pcontiguous=.01); filen="sinusStep10.jpg";mutationRate=2E-3;nbLocus=20;stepvalue=2;initial_locus_value=200;mutation_model="stepwise";mutation_rate=2E-2;nind="Ne";stat="F"
result=validation_with_coalescent(rasterStack=sinusStepStack,pK = pK, pR = pR,shapesK=shapesK,
                                  shapesR=shapesR, shapeDisp="contiguous_long_dist_mixt", 
                                  pDisp = c(pcontiguous=.05,plongdist=0.01), 
                                  filen="sinusStep10.jpg",mutationRate=2E-3,nbLocus=20, 
                                  stepvalue=2,initial_locus_value=200,mutation_model="stepwise",
                                  mutation_rate=2E-2,nind="Ne",stat="F")
rasterStack=sinusStepStack;pK = pK; pR = pR;shapesK=shapesK;shapesR=shapesR; shapeDisp="contiguous_long_dist_mixt"; pDisp = c(pcontiguous=.05,plongdist=0.01);filen="sinusStep10.jpg";mutationRate=2E-3;nbLocus=20; stepvalue=2;initial_locus_value=200;mutation_model="stepwise";mutation_rate=2E-2;nind="Ne";stat="F"
#transition = result$transitionmatrice

result=validation_with_coalescent(rasterStack=rasterStack,pK = pK, pR = pR,shapesK=shapesK,shapesR=shapesR, shapeDisp=shapeDisp, pDisp = pDisp, filen=NULL,mutationRate=1E-2,nbLocus=6, initial_locus_value=200)
result=validation_with_coalescent(rasterStack=rasterStack,pK = pK, pR = pR,shapesK=shapesK,shapesR=shapesR, shapeDisp=shapeDisp, pDisp = pDisp, filen="stepping10.jpg",mutationRate=1E-2,nbLocus=12, initial_locus_value=200)
plot_validation(result)

# small landscape
#
Data2 <- data.frame(BIO1=c(300,300,400,400),BIO12=c(2000,380,2900,2700)) 
N=150;r=10
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
#shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
#pR = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
shapesR=c(BIO1="constant",BIO12="constant")
pR = matrix(c(r,r),nrow=1,ncol=2,dimnames=list(c("Y"),c("BIO1","BIO12")))
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
smallStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=dim(Data2)[1]),xmn=0,xmx=dim(Data2)[1],ymn=0,ymx=1),"BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=dim(Data2)[1]),xmn=0,xmx=dim(Data2)[1],ymn=0,ymx=1)))
rasK=smallStack;values(rasK)= as.matrix(ReactNorm(values(smallStack),pK,shapesK)[,"Y"])
values(rasK)
result=validation_with_coalescent(rasterStack=smallStack,pK = pK, pR = pR,shapesK=shapesK,shapesR=shapesR, shapeDisp="contiguous8", pDisp = .05, filen=NULL,mutation_rate=1E-2,nbLocus=20, initial_locus_value=200,Option = "homogenous_1col",nind=10,stat="F")
plot_validation(result)
save(result,file="result.Rdata")
load("result.Rdata")
#
# test of Hitting times of Boley et al. 2012
#
P=matrix(0,nrow=6,ncol=6,dimnames=list(as.character(1:6),as.character(1:6)))
for (i in 1:5) {P[i,i+1]=1;P[6,1]=1;P[2,c(3,5)]=.5}
P;H=hitting_time_digraph(P);C=H+t(H)
C

#
# 2 dim random landscape
#
N=10
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
shapesR=c(BIO1="constant",BIO12="constant")
pR = matrix(c(r,r),nrow=1,ncol=2,dimnames=list(c("Y"),c("BIO1","BIO12")))
Data2 <- data.frame(BIO1=c(120,300,400,400),BIO12=c(380,1000,1300,2700)) 
Data2 <- data.frame(BIO1=sample(c(300,180,120,400),36,replace=TRUE),BIO12=sample(c(2000,1000,350,2900),36,replace=TRUE) )
#Data2 <- data.frame(BIO1=sample(300:101,36),BIO12=sample(31:230*10,36));dim(Data2);dimnames(Data2)
# on va fabriquer une raster stack avec deux couches correspondans aux variables environementales de la data frame
#rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=6,ncol=6),xmn=0,xmx=6,ymn=0,ymx=6),"BIO12"=raster(matrix(Data2$BIO12,nrow=6,ncol=6),xmn=0,xmx=6,ymn=0,ymx=6)))
#plot(rasterStack);writeRaster(rasterStack,"example.tif",format="GTiff",overwrite=TRUE)
rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=6,ncol=6),xmn=0,xmx=6,ymn=0,ymx=6),"BIO12"=raster(matrix(Data2$BIO12,nrow=6,ncol=6),xmn=0,xmx=6,ymn=0,ymx=6)))
rasK=rasterStack;values(rasK)= as.matrix(ReactNorm(values(rasterStack),pK,shapesK)[,"Y"])
values(rasK);plot(rasK)
writeRaster(rasterStack,"ExampleLAndscape6by6.tiff",foramt="Gtiff",overwrite=TRUE)
plot(r)
#rasterStack=rasterStack;pK = pK; pR = pR;shapesK=shapesK;shapesR=shapesR; shapeDisp="contiguous"; pDisp = .05; filen=NULL;mutationRate=1E-2;nbLocus=20; initial_locus_value=200;Option = "sample_1_col_diploid";nind=10;stat="F"

N=10
r=10
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
shapesR=c(BIO1="constant",BIO12="constant")
pR = matrix(c(r,r),nrow=1,ncol=2,dimnames=list(c("Y"),c("BIO1","BIO12")))
rasterStack <- stack("ExampleLAndscape6by6.tif")
names(rasterStack) <- c("BIO1","BIO12")
rasK=rasterStack;values(rasK)= as.matrix(ReactNorm(values(rasterStack),pK,shapesK)[,"Y"])
#rasterStack=rasterStack;pK = pK; pR = pR;shapesK=shapesK;shapesR=shapesR; shapeDisp="contiguous"; pDisp = .05; filen=NULL;mutationRate=1E-2;mutation_model="stepwise";nbLocus=20; initial_locus_value=200;Option = "homoNonEmtpyCells1col";Option = "homogenous_1col";nind="Ne2";stat="F";mutation_rate=1E-2;stepvalue=2
result=validation_with_coalescent(rasterStack=rasterStack,pK = pK, pR = pR,shapesK=shapesK,shapesR=shapesR, shapeDisp="contiguous", pDisp = .05, filen=NULL,mutation_rate=1E-2,nbLocus=40, initial_locus_value=200,Option = "homoNonEmtpyCells1col",nind="Ne2",stat="F",stepvalue=2)
plot_validation(result,whichCell="KmoreThan.5")
plot_validation(result,whichCell="KmoreThan.2",what=c(1,2,3,4,16,5,6,7,9),filen="example20150112.jpg")
coalTimesMatrixes <- checkTimeInterval(result$transitionmatrice,rasK)
plot_coal_time_depending_on_time_interval(coalTimesMatrixes,filen="coaltimes")


#
# Buseola fusca
#

environmentalData_ <- raster("AfrBfGAMActA.asc");environmentalData_[environmentalData_>200]=200
environmentalData_2 <- trim(environmentalData_, padding=(8*dim(environmentalData_)%/%8)[1:2]) 
environmentalData_2 <- aggregate(environmentalData_2,8)
genetData <- read.table("WBf16genelandcoord.txt")
genetData <- cbind(genetData,read.table("WBf16genelandgeno.txt"))
genetSP <- SpatialPoints(genetData[,c("x","y")])
bbox(genetSP)
environmentalData <- crop(environmentalData_2, bbox(genetSP)+res(environmentalData_2)*c(-2,-2,2,2))
genetData$Cell_numbers <- cellFromXY(environmentalData,genetData)
environmentalData[environmentalData==0] <- NA
ncellA(environmentalData)
any(is.na(extract(environmentalData,genetData[,c("x","y")])))
plot(environmentalData)
plot(genetSP,add=TRUE)
min(extract(environmentalData,genetData[,c("x","y")]))
genetData[which(extract(environmentalData,genetData[,c("x","y")])<0.01),]

genetData[is.na(genetData)] <- as.integer(1E9)
genetData=genetData[rowSums(genetData[,grep("ocus",colnames(genetData),value=TRUE)])<7E9,]
genetData[genetData==1E9]=NA
genetData <- TwoCols2OneCol(genetData)

# prior distributions
prior=list()
prior$K$distribution = "uniform"
prior$K$model = c(K="proportional")
prior$K$p = c(min=0.001,max=0.5)
prior$R$distribution = "fixed"
prior$R$model = c(K="constant")
prior$R$p = 20
prior$mutation_rate$model = "stepwise"
prior$mutation_rate$distribution = "loguniform"
prior$mutation_rate$p = c(min=1E-6,max=1E-2)
prior$dispersion$model="contiguous"
prior$dispersion$distribution="uniform"
prior$dispersion$p=c(min=0.001,max=0.5)

stepvalue = c(Locus1_65=3,Locus2_41=3,Locus3_117=2,Locus2_67=3,Locus1_70=3,Locus1_29=3,Locus3_21=2)

# stepping model
N=10;r=10
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
shapesR=c(BIO1="constant",BIO12="constant")
pR = matrix(c(r,r),nrow=1,ncol=2,dimnames=list(c("Y"),c("BIO1","BIO12")))
Data2 <- data.frame(BIO1=rep(120,10),BIO12=900,10);dim(Data2)
steppingStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=dim(Data2)[1]),xmn=0,xmx=dim(Data2)[1],ymn=0,ymx=1),"BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=dim(Data2)[1]),xmn=0,xmx=dim(Data2)[1],ymn=0,ymx=1)))
rasK=steppingStack;values(rasK)= as.matrix(ReactNorm(values(steppingStack),pK,shapesK)[,"Y"])
result=validation_with_coalescent(rasterStack=steppingStack,pK = pK, pR = pR,shapesK=shapesK,shapesR=shapesR, shapeDisp="contiguous", pDisp = .05, filen=NULL,mutation_rate=1E-2,nbLocus=40, initial_locus_value=200,Option = "homoNonEmtpyCells1col",nind="Ne2",stat="F",stepvalue=2)


# test simulation backward
# on va mettre des données environnementales de température et précipitations annuelles dans une data frame
#Data2 <- data.frame(BIO1=c(200,120,300,400),BIO12=c(1000,350,2000,2900)) 
Data2 <- data.frame(BIO1=c(300,200,200,400),BIO12=c(2000,380,320,2900)) 
Data2 <- data.frame(BIO1=sample(300:101,4),BIO12=sample(31:230*10,4));dim(Data2)
Data2 <- data.frame(BIO1=c(300,200,200,400),BIO12=c(2000,700,320,2900));dim(Data2)
# on va fabriquer une raster stack avec deux couches correspondans aux variables environementales de la data frame
rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1),"BIO12"=raster(matrix(Data2$BIO12,nrow=1,ncol=4),xmn=0,xmx=4,ymn=0,ymx=1)))
#rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=20,ncol=10),xmn=0,xmx=10,ymn=0,ymx=20),"BIO12"=raster(matrix(Data2$BIO12,nrow=20,ncol=10),xmn=0,xmx=10,ymn=0,ymx=20)))
plot(rasterStack)
# we test coalescence simulation fonction from invented coordinates values
# to have individuals sampled in areas where they are likely to be present
# we use the carrying capacity map
rasK=rasterStack;values(rasK)= as.matrix(ReactNorm(values(rasterStack),pK,shapesK)[,"Y"])
plot(rasK)
values(rasK)
geneticData = CreateGenetArray(rasK, 60,200,Option="full_1col_diploid")# the genetic data is meaning less (only one genotype), 
                                                              # but genetic data are not used. They will be modified according to the simulated coalescent
                                                              # we create as many individuals (lines) as K for each cell of the map 
                                                              # then we sample a few individuals for the coalescent
geneticData = CreateGenetArray(rasK, 60,200,Option="sample_1col_diploid",nind=8)
#geneticData = CreateGenetArray(rasK, 4,200,Option="sample_1col_diploid",nind=8)
geneticData$Cell_numbers
dim(geneticData)
geneticData$Cell_numbers
migrationMatrix(rasK,shapeDisp,pDisp)
coalgenet
system.time(coalescent_simulated <- simul_coalescent(geneticData,rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp,mutation_rate=1E-2,initial_genetic_value=200,mutation_model="stepwise",stepvalue=2,locusnames=NA))# we simulate a coalescent
coalescent_simulated$coalescent
apetree <- coalescent_2_newick(coalescent_simulated$coalescent)
# we plot the coalescent
plot_coalescent(coalescent_simulated,with_landscape=TRUE,rasK=rasK,legend_right_move=-.3)
plot_coalescent(coalescent_simulated,rasK=rasK,legend_right_move=-.1)
Qwithin_pair(coalescent_simulated$tip_genotype)
Qbetween(coalescent_simulated$tip_genotype)
aValue<-a_value_matrix(coalescent_simulated$tip_genotype,rasK)
aValueInd<-a_value_ind(coalescent_simulated$tip_genotype)
validation_with_coalescent(rasterStack=rasterStack,pK = pK, pR = pR,shapesK=shapesK,shapesR=shapeR, shapeDisp=shapeDisp, pDisp = pDisp, filen=NULL,mutationRate=1E-3,nbLocus=60, initial_locus_value=200)
validation_with_coalescent(rasterStack=rasterStack,pK = pK, pR = pR,shapesK=shapesK,shapesR=shapeR, shapeDisp="contiguous", pDisp = .05, filen=NULL,mutationRate=1E-3,nbLocus=60, initial_locus_value=200)

rasterStack <- stack(raster(matrix(120,ncol=10),xmn=0,xmx=10,ymn=0,ymx=1),raster(matrix(900,ncol=10),xmn=0,xmx=10,ymn=0,ymx=1))
names(rasterStack)<-c("BIO1","BIO12")
result=validation_with_coalescent(rasterStack=rasterStack,pK = pK, pR = pR,shapesK=shapesK,shapesR=shapesR, shapeDisp="contiguous", pDisp = .05, filen=NULL,mutationRate=1E-2,nbLocus=6, initial_locus_value=200)
plot_validation(result)
###

# check for larger data
Data2 <- data.frame(BIO1=sample(c(300,0,120,400),36,replace=TRUE),BIO12=sample(c(2000,0,350,2900),36,replace=TRUE) )
Data2 <- data.frame(BIO1=sample(300:101,36),BIO12=sample(31:230*10,36));dim(Data2);dimnames(Data2)
# on va fabriquer une raster stack avec deux couches correspondans aux variables environementales de la data frame
#rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=6,ncol=6),xmn=0,xmx=6,ymn=0,ymx=6),"BIO12"=raster(matrix(Data2$BIO12,nrow=6,ncol=6),xmn=0,xmx=6,ymn=0,ymx=6)))
#plot(rasterStack);writeRaster(rasterStack,"example.tif",format="GTiff",overwrite=TRUE)
rasterStack <- stack(("example.tif"));names(rasterStack)<-c("BIO1","BIO12")
Data2 <- values(rasterStack);dimnames(Data2)<- list(1:dim(Data2)[1],c("BIO1","BIO12"))
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
N=2;r=10
mutation_rate=1E-1
initial_genetic_value=200
mutation_model="stepwise"
stepvalue=2
locusnames=NA
pK = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
#shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
#pR = matrix(c(100,500,300,0,N,N,300,3000,2500,0,N,N),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
shapesR=c(BIO1="constant",BIO12="constant")
pR = matrix(c(r,r),nrow=1,ncol=2,dimnames=list(c("Y"),c("BIO1","BIO12")))
# we also have a shape of the dispersion function given here by the variable "shapeDisp"
# and dispersion parameters of the dispersion function
#pDisp = c(1/19,1)
pDisp = .5
#shapeDisp="fat_tail1"
shapeDisp="gaussian"
# we test coalescence simulation fonction from invented coordinates values
# to have individuals sampled in areas where they are likely to be present
# we use the carrying capacity map
rasK=rasterStack;values(rasK)= as.matrix(ReactNorm(values(rasterStack),pK,shapesK)[,"Y"])
plot(rasK)
# geneticData = CreateGenetArray(rasK, 20,200,Option="full_1col_diploid")
# the genetic data is meaning less (only one genotype), 
# but genetic data are not used. They will be modified according to the simulated coalescent
# we create as many individuals (lines) as K for each cell of the map 
# then we sample a few individuals for the coalescent
geneticData = CreateGenetArray(rasK, 60,200,Option="sample_1col_diploid",nind=100)
dim(geneticData)
geneticData$Cell_numbers
migrationMatrix(rasK,shapeDisp,pDisp)
# mutation_rate=1E-1,initial_genetic_value=200,mutation_model="stepwise",stepvalue=2,locusnames=NA
system.time(coalescent_simulated <- simul_coalescent(geneticData,rasterStack,pK,pR,shapesK,shapesR,shapeDisp,pDisp,mutation_rate=1E-2,initial_genetic_value=200,mutation_model="stepwise",stepvalue=2,locusnames=NA))# we simulate a coalescent
coalescent_simulated
apetree <- coalescent_2_newick(coalescent_simulated$coalescent)
# we plot the coalescent
plot_coalescent(coalescent_simulated,with_landscape=TRUE,rasK=rasK,legend_right_move=-.3,file="essai.pdf")
plot_coalescent(coalescent_simulated,rasK=rasK,legend_right_move=-.1)
aValueInd<-a_value_ind(coalescent_simulated$tip_genotype)

###
setwd(wd)
source("graphPOP_0.113.R")
#wd="/media/dupas/LACIE SHARE/IRD/ETUDIANTS/2014/Vasundra/R/" # fixe
Data2 <- na.omit(read.table("Data.csv")[c("BIO1","BIO12")])
Data2 <- data.frame(BIO1=rep((11:40)*10,10),BIO12=rep((0:9)*400,each=30))
rasterStack <- stack(list("BIO1"=raster(matrix(Data2$BIO1,nrow=30,ncol=10),xmn=0,xmx=10,ymn=0,ymx=30),"BIO12"=raster(matrix(Data2$BIO12,nrow=30,ncol=10),xmn=0,xmx=10,ymn=0,ymx=30)))
Data2 <- values(rasterStack)
pK = matrix(c(100,500,300,0,10,10,300,3000,2500,0,10,10),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
pR = matrix(c(100,500,300,0,10,10,300,3000,2500,0,10,10),nrow=6,ncol=2,dimnames=list(c("Xmin","Xmax","Xopt","Yxmin","Yxmax","Yopt"),c("BIO1","BIO12")))
shapesR=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
shapesK=c(BIO1="conquadraticskewed",BIO12="conquadraticskewed")
pDisp = c(3.6,1.6)
K = ReactNorm(values(rasterStack),pK,shapesK)[,"Y"]
r = ReactNorm(values(rasterStack),pR,shapesR)[,"Y"]
Data2 <- cbind(Data2,K)
rasK <- rasterStack;values(rasK) <- subset(Data2,select=K)
method="fat_tail1"
migrationM <- migrationMatrix(rasterStack,method, pDisp)
transitionmatrice = transitionMatrixBackward(r, K, migration= migrationM)
plot(raster(transitionmatrice))
geneticData = CreateGenetArray(rasK,  nb_locus=20, initial_locus_value=200,Option="full_1col_diploid")
geneticData = geneticData[sample(dim(geneticData)[1],100),];row.names(geneticData)=1:100
dim(geneticData)
coalescent_simulated <- simul_coalescent_simple(geneticData,transitionmatrice,dimGeneticData,rasK)
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
         r=list(shape="linear0",priorshape = "unif", pR=matrix(c(50,300,300,600,10,10,100,250,250,600,10,10),nrow=6,ncol=2,dimnames=list(c("rXY0Min","rXY0Max","rXmaxMin","rXmaxMax","rYXmaxMin","rYXmaxMax"),c("BIO1","BIO12"))),
                Tgen = list(shape="envelin",priorshape = "unif", pR=matrix(c(50,300,300,600,10,10),nrow=6,ncol=2,dimnames=list(c("rX00Min","rX00Max","rXmaxMin","rXmaxMax","rYXmaxMin","rYXmaxMax"),c("BIO1"))))))

matrix(c(100,500,300,0,10,10,100,500,300,0,10,10,100,500,60,30,300,3000,2500,0,10,10,300,3000,2500,0,10,10,300,3000,40,40),nrow=16,ncol=2,dimnames=list(c("KXmin","KXmax","KXopt","KYxmin","KYxmax","Yopt","rXmin","rXmax","rXopt","rYxmin","rYxmax","Yopt","GtimeXmin","GTimeXmax","GtimeYXmin","GTimeYXmax"),c("BIO1","BIO12")))
p_est = matrix(c(rep(TRUE,28),rep(FALSE,4)),nrow=16,ncol=2,dimnames=list(c("KXmin","KXmax","KXopt","KYxmin","KYxmax","Yopt","rXmin","rXmax","rXopt","rYxmin","rYxmax","Yopt","GtimeXmin","GTimeXmax","GtimeYXmin","GTimeYXmax"),c("BIO1","BIO12")))
pR = matrix(c(100,500,300,0,10,10,300,3000,2500,0,10,10),nrow=6,ncol=2,dimnames=list(c("rXmin","rXmax","rXopt","rYxmin","rYxmax","Yopt"),c("BIO1","BIO12")))
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
