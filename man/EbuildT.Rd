\name{EbuildT}
\alias{EbuildT}
\title{Builds expected values of sufficient statistics in multivariate
  normal regression}
\description{

  This does the E-step of an EM algorithm for multivariate regression.
  It is assumed that \code{X} is fully observed and \code{Y} may have
  missing data.  This returns an augmented matrix which contains the
  expected value of the means and the covariance matrix for the current
  set of parameters, \code{B}, \code{b0}, and \code{Syy.x}.


}
\usage{
MbuildT(X, Y, B, b0, Syy.x, w = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A data matrix with rows representing observations and columns
    variables.  These are independent variables in a regression, and
    must be fully observed.
  }
  \item{Y}{A data matrix with rows representing observations and columns
    variables.  These are dependent variables in a regression, and
    may contain missing values.  \code{X} and \code{Y} should have the
    same number of rows.
  }
  \item{B}{A \eqn{J \times K} matrix containing the non-constant
    coefficients of a multivariate regression of \code{Y} on \code{X}.
    Here \eqn{J} is the number of columns of \code{Y} and \eqn{K} is the
    number of columns of \code{X}.
  }
  \item{b0}{A vector of the constant terms from the multivariate
    regression of \code{Y} on \code{X}. This should have length
    \eqn{J}. 
  }
  \item{Syy.x}{A \eqn{J \times J} symmetric matrix giving the residual
    covariance of the multivariate regression of \code{Y} on \code{X}.
  }
  \item{w}{If supplied this should be a vector of weights equal to the
    number of rows of \code{X}.  The total sample size is considered to
    be the sum of the weights.
  }
}
\details{

  Following Little and Rubin (2002), the sufficient statistics for a
  multivariate normal distribution can be represented by a matrix:
  \deqn{ \bold{T} = \left [ \begin{array}{ccc}
    -1 & \mu_X & \mu_Y \\
    \mu_X & \Sigma_{XX} & \Sigma_{XY} \\
    \mu_Y & \Sigma_{YX} & \Sigma_{YY}
    \end{array} \right ] }.

  Sweeping (see \code{\link{matSweep}}) the rows and columns
  corresponding to the \eqn{\bold{X}} variables produces the
  multivariate regression of \eqn{\bold{Y}} on \eqn{\bold{X}}.
  
  \deqn{ SWP[X]\bold{T} = \left [ \begin{array}{ccc}
    * & * & b_{0.X} \\
    * & -\Sigma_{XX}^{-1} & B^T_{XY.X} \\
    b_{0.X} & B_{XY.X} & \Sigma_{YY.X}
    \end{array} \right ] },

  Here \eqn{B_{XY.X}} is the matrix of regression coefficients,
  \eqn{b_{0.X}} is the constants and \eqn{\Sigma_{YY.X}} is the residual
  covariances.

  If \code{Y} was fully observed, then \eqn{\bold{T}} could be easily
  calculated by first forming a matrix \code{XY1=cbind(1,X,Y)}, and then
  calculating \eqn{SWP[1](XY1)^T W (XY1)}, where \eqn{W} is a diagonal matrix
  with the weights, \code{w}, on the diagonal, and SWP[1] is the sweep
  operator applied to the constant row/column (assumed to be the first).
  Note that to get the parameter of the normal distribution, the part
  corresponding to the covariance matrix must be scaled by the sum of
  the weights.  

  If there are missing values, then calculating the expected value of
  \eqn{\bold{T}} requires two steps.  First, the missing values for
  \code{Y} must be imputed.  Second, an adjustment needs to be made for
  any term in the sum which involves the product of two imputed
  values (or the square of a single imputed value).

  In order to perform both the imputation and the adjustment, a 
  parameter estimate \eqn{bold{T}^{(i)}} is required.  In this
  function, this is supplied through the second parameterization,
  \eqn{SWP[X]\bold{T}}.  As \code{X} is fully observed the upper left
  \eqn{2 \times 2} submatrix can easily be calculated from the data.
  The remaining parts are passed in as arguments to the function.

  The principle output is \eqn{bold{T}^{(i+1)}}, that is one E-step
  of an EM-algorithm for the multivariate normal distribuiton.  Note
  that the imputed data matrix is also returned.

}
\value{
  A list with three components:
  \item{T}{A \eqn{(1+K+J) \times (1+K+J)} augmented covariance matrix
    estimate for \code{X} and \code{Y}.}
  \item{XY1}{The augmented data matrix with missing values of \code{Y}
    imputed  with regression imputations (using the old parameters).}
  \item{M}{A logical matrix showing were the missing values in \code{Y}
    are.}
}
\references{

    Little, R. J. A. and Rubin, D. B. (2002).  \emph{Statistical Analysis
    with Missing Data, Second Edition.}  Wiley.

}
\author{Russell Almond}
\note{
  This is built to facilitate \code{\link{BQreg}}.
}
\seealso{
  \code{\link{BQreg}}, \code{\link{matSweep}}.
}
\examples{

data(mvMux)

mvX <- mvMux[,1:2] # Complete Data
mvY <- mvMux[,3:5] # Data with missing values.

## Start with numbers based on complete data only.
mv.cccov <- cov.wt(na.omit(mvMux),method="ML")
mv.ccm <- apply(mvMux,2,mean,na.rm=TRUE)
mv.ccT <- rbind(c(-1/nrow(mvMux),mv.ccm),
                 cbind(mv.ccm,mv.cccov$cov))

## Use sweep to do the multivariate regression
mv.ccTswp12 <- matSweep(mv.ccT,2:3)
mv.BB <- mv.ccTswp12[4:6,2:3]
mv.B0 <- mv.ccTswp12[4:6,1]
mv.B <- cbind(mv.B0,mv.BB)

mv.Syy.x <- mv.ccTswp12[4:6,4:6]

Tstep <- EbuildT(mvX, mvY, mv.BB, mv.B0, mv.Syy.x, w=1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
