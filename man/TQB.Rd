\name{TQB}
\alias{TQB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produces a constrainted regression from an augmented covariance matrix.}
\description{

  A Q-constrained regression is a multivariate regression where the
  coefficients are constrated to be zero when the corresponding entries
  of the logical restriction matrix \code{Q} are false.  This is a
  utility function which produces the Q-constrained coefficients from
  the estimated augmented covariance matrix.


}
\usage{
TQB(T, Q)
}
\arguments{
  \item{T}{An augmented covariance matrix where the first row
    corresponds to the means, the next group of rows/columns the predictor
    variables and the final rows/columns the predictands.  It should be
    a square matrix of size \eqn{1+K+J} where \eqn{K} is the number of
    predictors and \eqn{J} is the number of predictands.
  }
  \item{Q}{A \eqn{J \times K} logical matrix where \eqn{K} is the number of
    predictors and \eqn{J} is the number of predictands.  This value
    should be \code{FALSE} if the correpsonding regression coefficient
    is constrained to zero, and \code{TRUE} otherwise.
  }
}
\details{

  The target is to produce a multivarite regression \eqn{\bold{Y} =
  \bold{X} \bold{B}+ \bold{b0} + \bold{E}_{Y.X}}, which is
  constrained by the matrix \eqn{\bold{Q}} so that if \eqn{q_{jk}=0}
  then \eqn{b_{jk} = 0}.

  It is possible to solve the regression using the sufficient statistics
  which   can be expressed as an augmented covariance matrix:
  \deqn{ \bold{T} = \left [ \begin{array}{ccc}
    -1 & \mu_X & \mu_Y \\
    \mu_X & \Sigma_{XX} & \Sigma_{XY} \\
    \mu_Y & \Sigma_{YX} & \Sigma_{YY}
    \end{array} \right ] }.
  The equations can then be solved with multiple applications of the
  sweep operator (\code{\link{matSweep}}).  

  If \eqn{\bold{X}} and \eqn{\bold{Y}} are fully observed, then
  \eqn{\bold{T}} can be easily esimated using
  \code{\link[SSX](cbind(X,Y))}.  If the variables are not fully
  observed then \eqn{\bold{T}} can be estimated using the EM-algorithm.
  (See \code{\link{EbuildT}} for the case where \eqn{\bold{Y}} has
  missing values.)

  This procedure also estimates the residual covariance matrix. Under
  the usual local independence assumption, the off-diagonal elements of
  this matrix should be zero.  Off diagonal elements are only calculated
  for pairs of Y variables which both have the same pattern in the
  Q-matrix.  Others are left at zero.

}
\value{
  A list with three components
  \item{B}{A numeric matrix of the same size as \code{Q} giving the regression
    coefficients.}
 \item{b0}{A numeric vector of length \code{nrow{Q}} which gives the
   constant terms.}
 \item{Syy.x}{A symmetric numeric matrix with the same number of rows as
   \code{Q} which provides the covariance of the residuals.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Russell Almond}
\note{

  This is really an interior function, and I probably shouldn't have
  exposed it, but doing so allows me to document that algorithm a little
  bit better.

}
\seealso{
  \code{\link{BQreg}}
}
\examples{

data(twoskill)  


## Estimated T matrix
That <- matSweep(SSX(as.matrix(twoskill.data0),1),1)
## Normalize covariance estimates
That[-1,-1] <- That[-1,-1]/nrow(twoskill.data0)

Best <- TQB(That,twoskill.Q)

## data0 has no noise, so should recover parameters exactly
stopifnot(all.equal(Best$B,twoskill.B,check.attributes=FALSE),
          all.equal(Best$b0,twoskill.b0,check.attributes=FALSE))

## data1 has a bit of residual standard deviation.
That1 <- matSweep(SSX(as.matrix(twoskill.data1),1),1)
## Normalize covariance estimates
That1[-1,-1] <- That1[-1,-1]/nrow(twoskill.data1)

Best1 <- TQB(That1,twoskill.Q)



}
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
\keyword{ regression }% __ONLY ONE__ keyword per line
