\name{Collisionijk}
\alias{Collisionijk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Collisionijk(Hitting_mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Hitting_mat}{
%%     ~~Describe \code{Hitting_mat} here~~
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
structure(function (Hitting_mat) 
{
    return(standardGeneric("Collisionijk"))
  }, generic = structure("Collisionijk", package = "GraphPOP"), package = "GraphPOP", group = list(), valueClass = character(0), signature = "Hitting_mat", default = `\001NULL\001`, skeleton = (function (Hitting_mat) 
stop("invalid call in method dispatch to 'Collisionijk' (no default method)", 
    domain = NA))(Hitting_mat), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
