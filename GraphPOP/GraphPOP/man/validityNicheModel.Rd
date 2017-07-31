\name{validityNicheModel}
\alias{validityNicheModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{validityNicheModel}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
validityNicheModel(object)
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
    if (class(object@variables) != "character")
        stop("error in NicheModel variables : variables just accept character!")
    if (!is.list(object@parameterList))
        stop("error in NicheModel parameterList : parameterList just accept list!")
    if (class(object@reactNorms) != "character")
        stop("error in NicheModel reactNorms : reactNorms just accept character!")
    if (FALSE \%in\% lapply(object@parameterList, is.numeric))
        stop("error in NicheModel parameter list : Parameter list just accept numeric!")
    if (length(object@variables) != length(object@reactNorms))
        stop("error in NicheModel : number of variables and number of reaction norms do not correspond")
    notMatching <- (unlist(lapply(1:length(object@parameterList),
        function(x) nbpar(object@reactNorms[x]) != length(object@parameterList[[x]]))))
    if (any(notMatching))
        stop(paste("error in NicheModel : number of paremeters and reactionNorm do not match for variable ",
            which(notMatching), ". ", sep = ""))
    TRUE
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
