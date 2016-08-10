\name{eggs}
\alias{eggs}
\docType{data}
\title{A kitched experiment involving the measurement of eggs}
\description{

  Dempster (1969) describes a kitchen experiment in which a dozen grade
  A eggs were measured with rough equipment.  The orginal measurement
  were the length, \eqn{L}, (in cm), width, \eqn{W}, (in cm) and volume,
  \eqn{V} (in cc)
}
\usage{data("eggs")}
\format{
  A matrix with 12 observations (row) and data on the following 3
  variables (columns):
  \describe{
    \item{logLength}{\eqn{log_{10} L}}
    \item{logWidth}{\eqn{log_{10} W}}
    \item{logVol}{\eqn{log_{10} (6/\pi) V}}
  }
}
\source{
  Dempsters, A. P. (1969).  \emph{Elements of Continuous Multivariate
  Analysis.}  Addison-Wesley.  Example 8.1, p. 151.

}
\examples{
data(eggs)

}
\keyword{datasets}
