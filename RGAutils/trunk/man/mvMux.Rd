\name{mvMux}
\alias{mvMux}
\docType{data}
\title{A permuted data set with missing values}
\description{

  A data set with missing values from Little and Rubin (2002), Table 7.4
  (p 154).  This version is permuted so that the complete data are the
  first two columns of the table.


}
\usage{data("mvMux")}
\format{
  A data frame with 13 observations on the following 5 variables.
  \describe{
    \item{\code{x3}}{a numeric vector (complete)}
    \item{\code{x5}}{a numeric vector (complete)}
    \item{\code{x1}}{a numeric vector (4 missing values)}
    \item{\code{x2}}{a numeric vector (4 missing values)}
    \item{\code{x4}}{a numeric vector (7 missing values)}
  }
}
\details{

  The data originally appeared in Draper and Smith (1981), and were used
  as example in Little and Rubin (2002, Table 7.4, p 154).  In the
  original example, \code{x5} was the dependent variable.

  As this was used as a test of a procedure for imputation with missing
  values in the columns have been permuted so that the two complete data
  columns appear first.  The column labels show the original order.
  
}
\source{
  These data in the original unpermuted form are available as
  \code{\link[mvnmle]{missvals}} in the \code{mvnmle} package.
}
\references{

  Draper, N.R. and Smith, H. (1981) Applied Regression Analysis. John
  Wiley and Sons. 

  Little, R. J. A. and Rubin, D. B. (2002).  \emph{Statistical Analysis
    with Missing Data, Second Edition.}  Wiley.


}
\examples{
data(mvMux)

}
\keyword{datasets}
