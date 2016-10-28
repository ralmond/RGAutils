\name{Categorical}
\alias{dcat}
\alias{rcat}
\alias{docat}
\alias{pocat}
\alias{qocat}
\alias{rocat}
\title{Categorical and Ordered Categorical Distributions.
}
\description{

    Probability and random number generator functions for categorical and ordered categorical data.  Data can be expressed as integers or factors.  For unordered categories, only \code{dcat} and \code{rcat} are supported.  For ordered categories, \code{pocat} and \code{qocat} are supported as well.

}
\usage{
dcat(x, prob, log = FALSE)
docat(x, prob, log = FALSE)
pocat(q, prob, log = FALSE)
qocat(p, prob, factor = TRUE, labels = ifelse(is.matrix(prob), colnames(prob), names(prob)))
rcat(n, prob, factor = TRUE, ordered = FALSE, labels = ifelse(is.matrix(prob), colnames(prob), names(prob)))
rocat(n, prob, factor = TRUE, labels = ifelse(is.matrix(prob), colnames(prob), names(prob)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{prob}{
%%     ~~Describe \code{prob} here~~
}
  \item{factor}{
%%     ~~Describe \code{factor} here~~
}
  \item{ordered}{
%%     ~~Describe \code{ordered} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
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
function (n, prob, factor = TRUE, ordered = FALSE, labels = ifelse(is.matrix(prob), 
    colnames(prob), names(prob))) 
{
    if (is.data.frame(prob)) 
        prob <- as.matrix(prob)
    if (is.matrix(prob)) {
        kmax <- ncol(prob)
        cp <- t(apply(prob, 1, cumsum))
        if (any(abs(cp[, kmax] - 1) > sqrt(.Machine$double.eps))) {
            stop("Probabilities don't sum to one.")
        }
    }
    else {
        kmax <- length(prob)
        cp <- cumsum(prob)
        if (abs(cp[kmax] - 1) > sqrt(.Machine$double.eps)) {
            stop("Probabilities don't sum to one.")
        }
    }
    u <- runif(n)
    if (is.matrix(cp)) {
        x <- apply(sweep(cp, 1, u, ">"), 1, sum) + 1
    }
    else {
        x <- apply(outer(u, cp, "<"), 1, sum) + 1
    }
    if (factor) {
        levels <- 1:max(x)
        if (is.null(labels)) {
            labels <- levels
        }
        x <- factor(x, levels, labels, ordered = ordered)
    }
    x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
