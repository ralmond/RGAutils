\name{orthoInvert}
\alias{orthoInvert}
\title{Inverts a matrix and finds the orthogolization transfomration}
\description{

  This takes a positive definite matrix \eqn{\bold{Q}} and finds a
  decomposition where \eqn{\bold{Q} = \bold{B}\bold{T}\bold{B}^T}, where
  \eqn{\bold{T}} is diagonal and \eqn{\bold{B}} is an upper triangular
  matrix.  At the same time if finds \eqn{\bold{A} = \bold{B}^{-1}} and
  the inverse of \eqn{\bold{Q}}.  

}
\usage{
orthoInvert(Q, tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{Q}{A symmetric positive definite matrix.  }
  \item{tol}{A threshold for finding zero roots.  If a diagonal element
    of \eqn{\bold{T}} is less than \code{tol} an error will be raised.
  }
}
\details{

  This follows from Dempster (1969).  Let \eqn{\bold{Q}} be the inner
  product matrix of a vector space (in statistical terms, this is better
  known as a covariance matrix).  Orthogonalization is a process that
  transverse the basis of that vector space into a series of orthoginal
  vectors (similar to principle components analysis).  Let \eqn{\bold{B}} be
  the upper triangular matrix which transforms the vectors.  Then 
  \eqn{\bold{Q} = \bold{B}\bold{T}\bold{B}^T}, where \eqn{\bold{T}} is
  a diagonal matrix.

  At the same time, note that \eqn{\bold{Q}^{-1} = \bold{A}^T
    \bold{T}^{-1} \bold{A}}, where \eqn{\bold{A}=\bold{B}^{-1}}.
  The \code{orthoInverse} function finds \eqn{\bold{A}}, \eqn{\bold{B}},
  \eqn{\bold{T}}, and \eqn{\bold{Q}^{-1}} by successive
  orthogonization.  


}
\value{
  A list with the following components:
  \item{T}{The diagonal matrix giving the lengths of the basis vectors
    in the new space.}
  \item{B}{The transformation matrix which translates the original basis
    to the orthoginal one.}
  \item{A}{The inverse transformation matrix which translates the
    orthoginal basis to the original one.}
  \item{Qinv}{The inverse of the argument.}
}
\references{
  Dempsters, A. P. (1969).  \emph{Elements of Continuous Multivariate
  Analysis.}  Addison-Wesley.
}
\author{Russell Almond}
\seealso{
  \code{\link{eigen}} and \code{\link{chol}} may be more useful for many
  purposes. 
}
\examples{
##Example is taken from Dempster (1969), pp 64--66.

## Initial data.
DQ <- matrix(c(19.1434,  9.0356,  9.7634, 3.2394,
                9.0356, 11.8658,  4.6232, 2.4746,
                9.7634,  4.6232, 12.2978, 3.8794,
                3.2394,  2.4746,  3.8794, 2.4604), 4,4)

## Expected results
DQT <- diag(c(19.1434, 7.60104, 7.318301, 1.117877))
DQB <- matrix(c(1,0.471995,0.510014,0.169218,
                0,1,0.0019629,0.1244067,
                0,0,1,0.304132,
                0,0,0,1),4,4)
DQA <- matrix(c(1,-0.471995,-0.509088,0.044331,
                0,1,-0.0019629,-0.123810,
                0,0,1,-0.304132,
                0,0,0,1),4,4)

oi <- orthoInvert(DQ)

stopifnot(all.equal(oi$T,DQT,tolerance=1e-4))
stopifnot(all.equal(oi$B,DQB,tolerance=1e-4))
stopifnot(all.equal(oi$A,DQA,tolerance=1e-4))
stopifnot(all.equal(oi$Qinv,solve(DQ)))

}
\keyword{ algebra }% use one of  RShowDoc("KEYWORDS")
\keyword{ array }% __ONLY ONE__ keyword per line
