\name{coalescent_2_newick}
\alias{coalescent_2_newick}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coalescent_2_newick(coalescent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coalescent}{
%%     ~~Describe \code{coalescent} here~~
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
structure(function (coalescent) 
{
    return(standardGeneric("coalescent_2_newick"))
  }, generic = structure("coalescent_2_newick", package = "GraphPOP"), package = "GraphPOP", group = list(), valueClass = character(0), signature = "coalescent", default = `\001NULL\001`, skeleton = (function (coalescent) 
stop("invalid call in method dispatch to 'coalescent_2_newick' (no default method)", 
    domain = NA))(coalescent), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
