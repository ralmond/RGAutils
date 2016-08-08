\name{SSX}
\alias{SSX}
\title{Compute weighted sums of variables, squares and cross products}
\description{

  This routine produces the sums of squares and cross products of an
  augmented data matrix.  The data matrix is first augmented with a
  column of \eqn{1}'s, so that in addition to the sums of squares and
  cross products, the sums of of the variables and the sum of the
  weights are calculated.

}
\usage{
SSX(X, weights = 1, side = "left")
}
\arguments{
  \item{X}{A data set expressed as a matrix (or an object which can be
    coerced into a matrix.  Rows represent individuals and columns
    represent variables.}

  \item{weights}{If supplied, this should should be a numeric vector of
  the same length as the number of rows in \code{X}.  This produces a
  weighted sum/sum of squares.}
  \item{side}{A character scalar which should have either the value
  \dQuote{left} or \dQuote{right} (abbreviations are allowed).  This
  controls on which side the of the matrix the column of \eqn{1}'s is
  added. }
}
\details{

  Dempster (1969) describes a trick where both the means and
  variances/covariances for a set of data can be calculated in a single
  matrix.  Let \eqn{\bold{X}} be a \eqn{n \times p} data matrix.  The
  augmented matrix, \eqn{\bold{X}_{(+)}}, is a \eqn{(n +1) \times p}
  matrix added by adding a column of \eqn{1}'s. If the augmented column
  is added on the right, then \eqn{\bold{X}_{(+)}^{T}\bold{X}_{(+)}} has
  the following form:
  \deqn{ \bold{Q}_{(+)} = \left [ \begin{array}{cc}
    \sum x_{i1}^2 & \sum x_{i1} x_{i2} & \ldots & \sum x_{i1} x_{ip} & \sum x_{i1}  \\
    \sum x_{i2} x_{i1} & \sum x_{i2}^2 & \ldots & \sum x_{i2} x_{ip} & \sum x_{i2} \\
    \vdots & \vdots & \ddots & \vdots & \vdots \\
    \sum x_{ip} x_{i1} & \sum x_{ip} x_{i2} & \ldots & \sum x_{ip}^2 & \sum x_{ip} \\
    \sum x_{i1} & \sum x_{i2} & \ldots & \sum x_{ip} & n
    \end{array} \right ]
  }
  
  Applying the sweep operator (\code{\link{matSweep}} has an interesting
  effect. It produces the matrix:

    \deqn{SWP[p+1]\bold{Q}_{(+)} = \left [ \begin{array}{cc}
    -\bold{T} & \bar \bold{X} \\
    \bar \bold{X} & -1/n
    \end{array} \right ] }

  where \eqn{\bold{T}} is \eqn{n-1} times the covariance matrix for
  \eqn{\bold{X}}.

  If weights are supplied, instead of calculating
  \eqn{\bold{X}_{(+)}^{T}\bold{X}_{(+)}}, \code{SSX} calculates
  \eqn{\bold{X}_{(+)}^{T}\bold{W}\bold{X}_{(+)}}, where \eqn{\bold{W}}
  is a diagonal matrix with the weights on the diagonal.  Thus it
  becomes a weighted sum of squares:
  \deqn{ \left [ \begin{array}{cc}
    \sum w_i x_{i1}^2 & \sum w_i x_{i1} x_{i2} & \ldots & \sum w_i x_{i1} x_{ip} & \sum w_i x_{i1}  \\
    \sum w_i x_{i2} x_{i1} & \sum w_i x_{i2}^2 & \ldots & \sum w_i x_{i2} x_{ip} & \sum w_i x_{i2} \\
    \vdots & \vdots & \ddots & \vdots & \vdots \\
    \sum w_i x_{ip} x_{i1} & \sum w_i x_{ip} x_{i2} & \ldots & \sum w_i x_{ip}^2 & \sum w_i x_{ip} \\
    \sum w_i x_{i1} & \sum w_i x_{i2} & \ldots & \sum w_i x_{ip} & \sum w_i
    \end{array} \right ]
  }

  While Dempster(1969) usually does the augmentation to the left, Little
  and Rubin (2002) does the augmentation to the left.  Other than the
  location of the augmented row and column in the final matrix, there is
  little difference.  The \code{side} argument controls on which side
  the augmentation is done.  The default is to follow the Little and
  Rubin style (on the left) rather than the on the right style shown in
  the equations above.

}
\value{
  A symmetric matrix of size \eqn{(p+1) \times (p+1)}, where \eqn{p} is
  the number of columns of the \code{X} argument.  Its contents are as
  described above.
}
\references{

  Dempsters, A. P. (1969).  \emph{Elements of Continuous Multivariate
  Analysis.}  Addison-Wesley.
  
  Little, R. J. A. and Rubin, D. B. (2002).  \emph{Statistical Analysis
    with Missing Data, Second Edition.}  Wiley.

}
\author{Russell Almond}

\seealso{
  \code{\link{matSweep}}, \code{\link[stats]{cov}},
  \code{\link[stats]{cov.wt}}, \code{\link[base]{crossprod}}
  
}
\examples{

data(eggs)
eggsQ <- SSX(eggs)
stopifnot(as.integer(eggsQ[1,1])==nrow(eggs))

## Sweep out constant row and column.
eggsT <- matSweep(eggsQ,1)
stopifnot(all(abs(eggsT[2:4,2:4]-(nrow(eggs)-1)*cov(eggs)) < .0001),
          all(abs(eggsT[1,2:4]-apply(eggs,2,mean)) < .0001),
          abs(eggsT[1,1]+1/nrow(eggs)) < .0001)

## Test of weighted version
wt <- c(rep(2/3,6), rep(1/3,6))
eggwQ <- SSX(eggs,wt,side="right")
stopifnot(abs(eggwQ[4,4]-sum(wt)) < .00001)

## Sweep
eggwT <- matSweep(eggwQ,4)
## Weighted mean.
eggwcov <- cov.wt(eggs,wt,method="ML")
stopifnot (all(abs(eggwT[1:3,1:3]/(sum(wt))- eggwcov$cov) < .0001),
           all(abs(eggwT[4,1:3]-eggwcov$center) < .0001))


}
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }% __ONLY ONE__ keyword per line
