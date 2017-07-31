\name{LandscapeHistory}
\alias{LandscapeHistory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{LandscapeHistory}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LandscapeHistory(Landscapelist = listOfLandscape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Landscapelist}{
%%     ~~Describe \code{Landscapelist} here~~
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
function (Landscapelist = listOfLandscape)
{
    li <- unlist(lapply(1:length(Landscapelist), function(x) Landscapelist[[x]]["period"][1]))
    o <- order(c(li))
    lo <- unlist(lapply(1:length(o), function(x) lapply(1:length(o),
        function(y) {
            if (o[y] == x)
                Landscapelist[y]
        })))
    new("LandscapeHistory", lo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
