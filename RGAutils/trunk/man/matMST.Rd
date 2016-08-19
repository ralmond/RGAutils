\name{matMST}
\alias{matMST}
\title{Finds a partial orthonormal basis for a inner product matrix.}
\description{

  The multistandarize procedure is a method for finding an orthonormal
  basis for a given inner product (covariance) matrix.  In a series of
  steps it transforms an inner product matrix into the identity matrix,
  producing the matrix that would transform the vector space onto the
  orthonormal basis in the process.

}
\usage{
matMST(Q, C = diag(nrow(Q)), do = 1L:nrow(Q), done = integer())
}
\arguments{
  \item{Q}{A symmetric positive definite matrix.}
  \item{C}{Normally this is the identity matrix with the same number of
    rows and columns as \code{Q}, but other matrixes can be used as
    well.  See details.}
  \item{do}{A vector of integers indicating the row/columns of \code{Q}
    that are to be standardized.}
  \item{done}{A vector of integers indicating the row/columns of
    \code{Q} which have already been standardized.}
}
\details{

  This operator was originally proposed by Beaton (1964) and is
  described in Dempster (1969, pp 70-74).  Standardizing a variable
  basically ensures that its variance (i.e., inner product) is 1.
  Orthogonolization produces a basis that is orthogonol (i.e., the
  variables are uncorrelated).

  The procedure is iterative and gradually replaced \code{Q} with the
  identity matrix and \code{C} with a matrix which can be used to
  transform variables onto the orthonormal basis.  The \code{do}
  argument allows finer control over which variables are subjected to
  the orthonormalizaiton process.  In particular, if only the first
  \eqn{k} variables are transformed, then the upper \eqn{k \times k}
  submatrix of \code{Q} will be the identity matrix.  Also, the
  submatrix with the remaining rows and columns of \code{C} will still
  be the identity matrix (if \code{C} was the identity matrix to start).

  The result of applying \code{matMST} to a subset of rows can be used
  as the input to a new call to \code{matMST} to complete the process.
  In that case, the \code{done} argument should be used to indicate
  which variables have already been standardized.

  Finally, \code{C} does not need to be the identity matrix.  In
  particular, if \code{Q} is already from a transformed space, where
  \eqn{\bold{K}} is the transformation matrix which produced the
  variables used in \eqn{\bold{Q}} then inputing \eqn{\bold{K}} as the
  arugment \code{C} will produce \eqn{\bold{C}\bold{K}} as the output
  (where \eqn{\bold{C}} was the transformation produced with the
  identity matrix as input).

}
\value{
  A list with the following components.
  \item{Q}{The matrix \eqn{\bold{Q}_{do,done}} where the variables
    corresponding to \code{c(do,done)} have been normalized.}
  \item{C}{The transformation which will transform theoriginal variable
    space to the partially orthnomalized basis.}
  \item{done}{An integer vector indicating which variables have been
    orthonormalized.}
}
\references{

  Beaton, A. E. (1964).  The Use of Special Matrix Operators in
  Statistical Calculus.  Educational Testing Service Research Report
  RB-64-51. 

  Dempsters, A. P. (1969).  \emph{Elements of Continuous Multivariate
  Analysis.}  Addison-Wesley.

}
\author{Russell Almond}
\note{

  The orthogonolization is not unique; hence I haven't seen this used
  much in practice.  In particular, the output is order dependent, as
  the order will influence the orthogonolization.

  In particular, note that the standard applicaiton of standardizing on
  all of the columns is basically the inverse of the Cholesky
  decomposition.  In particular, \code{\link{solve}(\link{chol}(Q))} is
  probababy faster and more numerically stable.
  
}
\seealso{
  \code{\link{orthoInvert}}, \code{\link{chol}}
}
\examples{
## Test of Multistandardize, follows Dempster (1969, p. 71-3)
DQ <- matrix(c(19.1434,  9.0356,  9.7634, 3.2394,
                9.0356, 11.8658,  4.6232, 2.4746,
                9.7634,  4.6232, 12.2978, 3.8794,
                3.2394,  2.4746,  3.8794, 2.4604), 4,4)

DQmst <- matMST(DQ)

stopifnot(all.equal(DQmst$Q,diag(4)))
stopifnot(all.equal(DQmst$C,t(solve(chol(DQ)))))
stopifnot(all.equal(DQmst$done,1:4))

}
\keyword{ algebra }% use one of  RShowDoc("KEYWORDS")
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
