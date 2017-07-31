\name{validityMigrationModel}
\alias{validityMigrationModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{validityMigrationModel}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
validityMigrationModel(object)
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
function (object)
{
    if (!is.character(object@shapeDisp))
        stop("error in  MigrationModel shapeDisp : ShapeDisp just accept character!")
    if (length(object@shapeDisp) != 1)
        stop("error in  MigrationModel shapeDisp : ShapeDisp is  not valid because it contains more or less than one shape!")
    if (!object@shapeDisp \%in\% listOfMigrationShape)
        stop("error in  MigrationModel shapeDisp : the given shape not exist for MigrationModel !")
    if (FALSE \%in\% lapply(object@pDisp, is.numeric))
        stop("error in MigrationModel pDisp : pDisp just accept numeric!")
    if (nbpar(object@shapeDisp) != length(object@pDisp))
        stop("error in MigrationModel : number of paremeters and shapeDisp do not match for variable")
    TRUE
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
