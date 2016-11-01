\name{median.ordered}
\alias{median.ordered}
\title{The median operation for ordered categorical variables.
}
\description{

  This method extends the definition to ordered factors.
}
\usage{
\S3method{median}{ordered}(x, na.rm = FALSE)
}
\arguments{
  \item{x}{A vector of ordered categorical values.}
  \item{na.rm}{A logical value.  If true, then \code{NA}s are removed
    prior to calculating the median.}
}
\details{

  Although medians are not well defined for nominal scale variables,
  they are well defined for ordinal scale variables (i.e.,
  \code{\link[base]{ordered}} factors).

  If the vector length is odd, then \code{\link[base]{round}()} is used
  to determine whether the upper or lower endpoint of the interval is used.

}
\value{
  A ordered factor of length one.
}
\author{Russell Almond}
\seealso{
  \code{\link{qocat}}
}
\examples{
x <- ordered(sample.int(3,10,replace=TRUE),levels=1:3,labels=c("L","M","H"))
median(x)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
