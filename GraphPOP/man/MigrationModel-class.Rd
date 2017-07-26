\name{MigrationModel-class}
\Rdversion{1.1}
\docType{class}
\alias{MigrationModel-class}
\alias{[,MigrationModel,character,missing,ANY-method}
\alias{migrationMatrix,Landscape,MigrationModel-method}

\title{Class \code{"MigrationModel"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MigrationModel", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{shapeDisp}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{pDisp}:}{Object of class \code{"ANY"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "MigrationModel", i = "character", j = "missing", drop = "ANY")}: ... }
    \item{migrationMatrix}{\code{signature(object = "Landscape", model = "MigrationModel")}: ... }
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
showClass("MigrationModel")
}
\keyword{classes}
