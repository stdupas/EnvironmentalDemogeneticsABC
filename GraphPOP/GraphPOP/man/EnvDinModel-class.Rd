\name{EnvDinModel-class}
\Rdversion{1.1}
\docType{class}
\alias{EnvDinModel-class}
\alias{[,EnvDinModel,character,missing,ANY-method}
\alias{createDemographic,Landscape,EnvDinModel-method}
\alias{runEnvDinModel,Landscape,EnvDinModel-method}

\title{Class \code{"EnvDinModel"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EnvDinModel", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{K}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{R}:}{Object of class \code{"ANY"} ~~ }
    \item{\code{migration}:}{Object of class \code{"ANY"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "EnvDinModel", i = "character", j = "missing", drop = "ANY")}: ... }
    \item{createDemographic}{\code{signature(object = "Landscape", model = "EnvDinModel")}: ... }
    \item{runEnvDinModel}{\code{signature(object = "Landscape", model = "EnvDinModel")}: ... }
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
showClass("EnvDinModel")
}
\keyword{classes}
