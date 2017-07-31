\name{ordinary_laplacian}
\alias{ordinary_laplacian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ordinary_laplacian}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ordinary_laplacian(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
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
structure(function (object)
{
    return(standardGeneric("ordinary_laplacian"))
  }, generic = structure("ordinary_laplacian", package = "GraphPOP"), package = "GraphPOP", group = list(), valueClass = character(0), signature = "object", default = `\001NULL\001`, skeleton = (function (object)
stop("invalid call in method dispatch to 'ordinary_laplacian' (no default method)",
    domain = NA))(object), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
