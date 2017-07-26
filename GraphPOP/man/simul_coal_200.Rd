\name{simul_coal_200}
\alias{simul_coal_200}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simul_coal_200(demographic, printCoal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{demographic}{
%%     ~~Describe \code{demographic} here~~
}
  \item{printCoal}{
%%     ~~Describe \code{printCoal} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (demographic, printCoal) 
{
    return(standardGeneric("simul_coal_200"))
  }, generic = structure("simul_coal_200", package = "GraphPOP"), package = "GraphPOP", group = list(), valueClass = character(0), signature = c("demographic", 
"printCoal"), default = `\001NULL\001`, skeleton = (function (demographic, 
    printCoal) 
stop("invalid call in method dispatch to 'simul_coal_200' (no default method)", 
    domain = NA))(demographic, printCoal), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
