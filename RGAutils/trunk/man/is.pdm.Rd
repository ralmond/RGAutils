\name{is.pdm}
\alias{is.pdm}
\title{Checks for postivie definite matrixes}
\description{
  Tests to see if a matrix is symmetric and positive definite.
}
\usage{
is.pdm(m, tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{m}{A matrix.  (Function returns \code{FALSE} if this argument is
    not a square matrix.)}
  \item{tol}{A tolerance value for the symmetry test and in the check for 0
    eigenvalues. }
}
\details{
  An object \code{m} is a positive definite matrix if it satisfies the
  following tests in order:
  \enumerate{
    \item{It is a numeric matrix with no missing values.}
    \item{It is a square matrix.}
    \item{It is symmetric:  when subtracted from its transpose, the
      difference is not greater than \code{tol} in absolute value.}
    \item{There are no zero eigenvalues (eigenvalues with values less
      than \code{tol}).}
  }

  Positive definite matrixes are important because they can be used as
  covariances matrixes, among other uses.
}
\value{
  A logical value indicating whether or not the first argument is a
  positive definite matrix.
}
\author{Russell Almond}
\note{
  The package \code{matrixcalc} has a function called
  \code{is.positive.definite} which is similar.  There are a couple of
  differences.  First, this version uses a more numerically stable check
  for symmetry.  Second, under some conditions (e.g., \code{NA}s in the
  matrix) this version returns \code{FALSE} instead raising an error.
}
\seealso{
  \code{\link[base]{eigen}}
}
\examples{

## Diagonal matrixes are always positive definite
stopifnot(is.pdm(diag(3))

stopifnot(!is.pdm(matrix(1,3,3)))


}
\keyword{ algebra }% use one of  RShowDoc("KEYWORDS")
\keyword{ array }% use one of  RShowDoc("KEYWORDS")

