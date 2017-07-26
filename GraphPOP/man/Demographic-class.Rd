\name{Demographic-class}
\Rdversion{1.1}
\docType{class}
\alias{Demographic-class}
\alias{Demographic}
\alias{[,Demographic,character,missing,ANY-method}
\alias{compare,Demographic,Landscape,logical-method}
\alias{nCellA,Demographic-method}
\alias{simul_coal_200,Demographic,logical-method}
\alias{simul_coalescent,Demographic,logical-method}

\title{Class \code{"Demographic"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{Demographic(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{K}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{R}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{TransiBackw}:}{Object of class \code{"TransitionBackward"} ~~ }
    \item{\code{TransiForw}:}{Object of class \code{"TransitionForward"} ~~ }
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
Class \code{"\linkS4class{Landscape}"}, directly.
Class \code{"\linkS4class{RasterStack}"}, by class "Landscape", distance 2.
Class \code{"\linkS4class{Raster}"}, by class "Landscape", distance 3.
Class \code{"\linkS4class{RasterStackBrick}"}, by class "Landscape", distance 3.
Class \code{"\linkS4class{BasicRaster}"}, by class "Landscape", distance 4.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "Demographic", i = "character", j = "missing", drop = "ANY")}: ... }
    \item{compare}{\code{signature(demographic = "Demographic", popSize = "Landscape", printCoal = "logical")}: ... }
    \item{nCellA}{\code{signature(object = "Demographic")}: ... }
    \item{simul_coal_200}{\code{signature(demographic = "Demographic", printCoal = "logical")}: ... }
    \item{simul_coalescent}{\code{signature(demographic = "Demographic", printCoal = "logical")}: ... }
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
showClass("Demographic")
}
\keyword{classes}
\usage{
Demographic(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
