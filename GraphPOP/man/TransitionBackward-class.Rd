\name{TransitionBackward-class}
\Rdversion{1.1}
\docType{class}
\alias{TransitionBackward-class}
\alias{commute_time_digraph,TransitionBackward-method}
\alias{commute_time_undigraph,TransitionBackward-method}
\alias{hitting_time_digraph,TransitionBackward-method}
\alias{laplaceMatrix,TransitionBackward-method}
\alias{linearizedFstDigraph,TransitionBackward,Landscape-method}
\alias{linearizedFstUndigraph,TransitionBackward,Landscape-method}
\alias{ordinary_laplacian,TransitionBackward-method}

\title{Class \code{"TransitionBackward"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TransitionBackward", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{matrix}"}, from data part.
Class \code{"\linkS4class{array}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{mMatrix}"}, by class "matrix", distance 2.
Class \code{"\linkS4class{structure}"}, by class "matrix", distance 3.
Class \code{"\linkS4class{vector}"}, by class "matrix", distance 4, with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{commute_time_digraph}{\code{signature(object = "TransitionBackward")}: ... }
    \item{commute_time_undigraph}{\code{signature(object = "TransitionBackward")}: ... }
    \item{hitting_time_digraph}{\code{signature(object = "TransitionBackward")}: ... }
    \item{laplaceMatrix}{\code{signature(object = "TransitionBackward")}: ... }
    \item{linearizedFstDigraph}{\code{signature(transition = "TransitionBackward", popSize = "Landscape")}: ... }
    \item{linearizedFstUndigraph}{\code{signature(transition = "TransitionBackward", popSize = "Landscape")}: ... }
    \item{ordinary_laplacian}{\code{signature(object = "TransitionBackward")}: ... }
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
showClass("TransitionBackward")
}
\keyword{classes}
