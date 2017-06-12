
genetic <- setClass("genetic",
                    slots = c(ploidy="integer",ploidyByrow="logical"),
                    contains="data.frame",
                    prototype = prototype(data.frame(locus1.1=c(200,202),locus1.2=c(204,200)),ploidy=as.integer(2), ploidyByrow=FALSE),
                    validity = function(object){
                      if (all(grepl("ocus",names(object)))) TRUE else stop("col names of genetic data.frame do not contain 'ocus'")
                      if ((object@ploidy==2)&(object@ploidyByrow==FALSE)) {
                        if (length(grep("\\.1",names(object)))==0|length(grep("\\.2",names(object)))==0) {
                          if ((grep("\\.1",names(object))%%2!=1)|(grep("\\.2",names(object))%%2!=0)){
                            stop("Columns of diploid by row FALSE data frame have to be named as follows: 'c('.1','.2','.1','.2')'")
                          }
                        }
                      }
                    }
                    )

genetic <- function(df=data.frame(locus1.1=c(200,202),locus1.2=c(204,200)),ploidy=NA, ploidyByrow=NA){
  if (is.na(ploidy)) { 
    if (any(grep("\\.2", colnames(df)))) 
      {
      P <- (unlist(strsplit(grep("\\.",colnames(df),value=TRUE),"\\.")))
      ploidy <- max(suppressWarnings(as.integer(P))[!is.na(suppressWarnings(as.integer(P)))])
    }
    if (any(grep("\\.2", rownames(df)))) 
    {
      P <- (unlist(strsplit(grep("\\.",rownames(df),value=TRUE),"\\.")))
      ploidy <- max(suppressWarnings(as.integer(P))[!is.na(suppressWarnings(as.integer(P)))])
    }
  }
  if (is.na(ploidyByrow)) ploidyByrow = !(any(grep(paste("\\.",ploidy,sep=""), colnames(df))))
  new("genetic",df,ploidy=ploidy,ploidyByrow=ploidyByrow)
}

spatialGenetic <- setClass("spatialGenetic",
                           slots = c(x="numeric", y="numeric",Cell_numbers="numeric"),
                           contains = "genetic",
                           prototype = prototype(genetic(),x=c(1,2),y=c(1,1),Cell_numbers=c(1,2)),
                           validity = function(object){
                             if (length(object@x)!=length(object@y)) stop("slots x and y do not have the same length")
                             if (length(object@x)!=nrow(object)) stop("slots x and genetic do not have the same number of individuals")
                           }
                           )

spatialGenetic <- function(df=NA,x=NA,y=NA,Cell_numbers=NA)
{
  if(is.na(df)) df=cbind(data.frame(genetic()),x=c(1,2),y=c(1,1),Cell_numbers=c(1,2))
  if (!(all(c("x","y","Cell_numbers")%in%colnames(df)))){
      if (is.na(Cell_numbers)){
        if ("Cell_numbers"%in%colnames(df)) {
          Cell_numbers=df$Cell_numbers
          df <- df[,-which(colnames(df)=="Cell_numbers")]
        }
      }
      if (is.na(x)){
        if ("x"%in%colnames(df)) {
          x=df$x
          df <- df[,-which(colnames(df)=="x")]
        }
      }
      if (is.na(y)){
        if ("x"%in%colnames(df)) {
          y=df$y
          df <- df[,-which(colnames(df)=="y")]
        }
      }
    }
new("spatialGenetic",genetic(df[,-which(colnames(df)%in%c("x","y","Cell_numbers"))]),x=df$x,y=df$y,Cell_numbers=df$Cell_numbers)
}
  
