\name{Landscape-class}
\Rdversion{1.1}
\docType{class}
\alias{Landscape-class}
\alias{[,Landscape,character,missing,ANY-method}
\alias{compare,Demographic,Landscape,logical-method}
\alias{createDemographic,Landscape,EnvDinModel-method}
\alias{linearizedFstDigraph,TransitionBackward,Landscape-method}
\alias{linearizedFstUndigraph,TransitionBackward,Landscape-method}
\alias{migrationMatrix,Landscape,MigrationModel-method}
\alias{runEnvDinModel,Landscape,EnvDinModel-method}
\alias{runNicheModel,Landscape,NicheModel-method}
\alias{transitionMatrixBackward,Landscape,list-method}
\alias{valuesA,Landscape-method}

\title{Class \code{"Landscape"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Landscape", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{period}:}{Object of class \code{"Date"} ~~ }
    \item{\code{vars}:}{Object of class \code{"character"} ~~ }
    \item{\code{distanceMatrix}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{filename}:}{Object of class \code{"character"} ~~ }
    \item{\code{layers}:}{Object of class \code{"list"} ~~ }
    \item{\code{title}:}{Object of class \code{"character"} ~~ }
    \item{\code{extent}:}{Object of class \code{"Extent"} ~~ }
    \item{\code{rotated}:}{Object of class \code{"logical"} ~~ }
    \item{\code{rotation}:}{Object of class \code{".Rotation"} ~~ }
    \item{\code{ncols}:}{Object of class \code{"integer"} ~~ }
    \item{\code{nrows}:}{Object of class \code{"integer"} ~~ }
    \item{\code{crs}:}{Object of class \code{"CRS"} ~~ }
    \item{\code{history}:}{Object of class \code{"list"} ~~ }
    \item{\code{z}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{RasterStack}"}, directly.
Class \code{"\linkS4class{Raster}"}, by class "RasterStack", distance 2.
Class \code{"\linkS4class{RasterStackBrick}"}, by class "RasterStack", distance 2.
Class \code{"\linkS4class{BasicRaster}"}, by class "RasterStack", distance 3.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "Landscape", i = "character", j = "missing", drop = "ANY")}: ... }
    \item{compare}{\code{signature(demographic = "Demographic", popSize = "Landscape", printCoal = "logical")}: ... }
    \item{createDemographic}{\code{signature(object = "Landscape", model = "EnvDinModel")}: ... }
    \item{linearizedFstDigraph}{\code{signature(transition = "TransitionBackward", popSize = "Landscape")}: ... }
    \item{linearizedFstUndigraph}{\code{signature(transition = "TransitionBackward", popSize = "Landscape")}: ... }
    \item{migrationMatrix}{\code{signature(object = "Landscape", model = "MigrationModel")}: ... }
    \item{runEnvDinModel}{\code{signature(object = "Landscape", model = "EnvDinModel")}: ... }
    \item{runNicheModel}{\code{signature(object = "Landscape", model = "NicheModel")}: ... }
    \item{transitionMatrixBackward}{\code{signature(object = "Landscape", model = "list")}: ... }
    \item{valuesA}{\code{signature(object = "Landscape")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("Landscape")
}
\keyword{classes}
