\name{simul_multi_coal}
\alias{simul_multi_coal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{simul_multi_coal}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simul_multi_coal(demographic, printCoal,iteration)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{demographic}{
%%     ~~Describe \code{demographic} here~~
}
  \item{printCoal}{
%%     ~~Describe \code{printCoal} here~~
}
  \item{iteration}{
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
structure(function (demographic, printCoal,iteration)
{
    return(standardGeneric("simul_multi_coal"))
  }, generic = structure("simul_multi_coal", package = "GraphPOP"), package = "GraphPOP", group = list(), valueClass = character(0), signature = c("demographic",
"printCoal","numeric"), default = `\001NULL\001`, skeleton = (function (demographic,
    printCoal,iteration)
stop("invalid call in method dispatch to 'simul_multi_coal' (no default method)",
    domain = NA))(demographic, printCoal, iteration), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
