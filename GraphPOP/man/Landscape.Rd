\name{Landscape}
\alias{Landscape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Landscape}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Landscape(rasterstack = rasterstack, period = dateVector, vars = charVector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rasterstack}{
%%     ~~Describe \code{rasterstack} here~~
}
  \item{period}{
%%     ~~Describe \code{period} here~~
}
  \item{vars}{
%%     ~~Describe \code{vars} here~~
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
function (rasterstack = rasterstack, period = dateVector, vars = charVector)
{
    if (class(rasterstack) != "RasterStack")
        stop("error in Landscape rasterstack : rasterstack just accept RasterStack!")
    if (length(rasterstack@layers) == 0)
        stop("rasterstack values is null!")
    if (class(period) != "Date")
        stop("error in Landscape period : period just accept Date!")
    if (!is.character(vars))
        stop("error in Landscape vars : vars just accept character!")
    if (length(period) == 1)
        period <- c(period, period)
    if (length(names(rasterstack)) != length(vars))
        stop("error when creating landscape : the number of layers in the rasterstack differs from the number of variables in var")
    names(rasterstack) <- vars
    b <- xyFromCellA(rasterstack)
    mat = sapply(1:nrow(b), function(l1) {
        sapply(1:nrow(b), function(l2) {
            if (l1 == l2)
                sqrt(res(rasterstack)[1]^2 + res(rasterstack)[2]^2)/3
            else {
                sqrt((b[l1, 1] - b[l2, 1])^2 + (b[l1, 2] - b[l2,
                  2])^2)
            }
        })
    })
    new("Landscape", rasterstack, period = period, vars = vars,
        distanceMatrix = mat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
