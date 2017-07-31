\name{nbpar}
\alias{nbpar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{nbpar}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nbpar(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x)
{
    unlist(lapply(x, function(x) switch(x[], constant = 1, proportional = 1,
        enveloppe = 2, envelin = 2, quadratic = 4, fat_tail1 = 2,
        gaussian = 1, exponential = 1, contiguous = 1, contiguous8 = 1,
        island = 1, fat_tail2 = 2, contiguous_long_dist_mixt = 2,
        gaussian_long_dist_mixt = 2)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
