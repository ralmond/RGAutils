\name{BQreg}
\alias{BQreg}
\title{Multivariate regression restricted by a Q-matrix}
\description{

  This performs a multivariate regression where the coefficients are
  restricted by a \eqn{Q}-matrix.  In particular, the elements of the
  coefficient matrix corresponding to false entries in the
  \eqn{Q}-matrix are forced to be zero.  The motivation is a linear
  psychometric model, where the \code{X} variables represent latent
  abilities, the \code{Y} variables represent observables and only
  certain abilities are relevant for certain observables.

}
\usage{
BQreg(X, Y, Q = matrix(TRUE, ncol(Y), ncol(X)), weights = 1, maxstep = 25L, tolerance = 10 * sqrt(.Machine$double.eps), trace = FALSE)
}
\arguments{
  \item{X}{A \eqn{N \times K} matrix of predictor variables.  This is
    assumed to be fully observed.  As a special case, if \code{X} is a
    data frame with numeric entries, it is coverted to a matrix.
  }
  \item{Y}{A \eqn{N \times J} matrix of observable outcome variables.
    Missing data are allows in this part of the equation. As a special
    case, if \code{Y} is a data frame with numeric entries, it is
    coverted to a matrix. 
  }
  \item{Q}{If supplied, this should be a \eqn{J \times K} logical matrix
    where a \eqn{q_{jk}} is \code{TRUE} exactly when predictor \eqn{k} is
    relevant for observation \eqn{j}.
  }
  \item{weights}{This should either be a numeric scalar, or a vector of
    length \eqn{N} providing the weight.  The sum of the weights is
    taken as the sample size.
  }
  \item{maxstep}{An integer scalar.  If the EM-algorithm is used because
    of missing data, the maximum number of EM-cycles to use.
  }
  \item{tolerance}{A numeric scalar.  If the EM-algorithm is used
    because of missing data,  the threshold for convergence.
  }
  \item{trace}{A logical scalar.  If true, tracing output is printed
    during the EM algorithm.
  }
}
\details{

  The target is to produce a multivarite regression \eqn{\bold{Y} =
  \bold{X} \bold{B}+ \bold{b0} + \bold{E}_{Y.X}}, which is
  constrained by the matrix \eqn{\bold{Q}} so that if \eqn{q_{jk}=0}
  then \eqn{b_{jk} = 0}.  This is motivated by cognitively diagnostic
  psychometric models where the \eqn{Q}-matrix is often used to express
  the relationship between a collection of (usually latent) proficiency
  variables and a collection of observable outcome variables.  (See
  \code{\link{TQB}} for this part of the function.)

  The motivation is to use this constrained regression as part of a
  latent variable model where the \code{X} variables are estimated
  through importance sampling.  Therefore, there may also be weights
  (supplied through the variable \code{weights}).  In this case, a
  weighted least squares is performed.

  An additional issue is that one or more of the \code{Y} values may be
  missing.  (As the \code{X} values are assumed to be imputed, they are
  considered fully observed.)  In this case the EM-algorithm is applied
  to find an approximate maximum likelihood estimate.  Two assumptions
  are necessary for this algorithm to work:  (1) the \code{Y} values are
  missing at random and (2) the \code{Y} values are approximately
  normally distributed.  Note that the matrix sweep operator
  (\code{\link{matSweep}}) is used for as much of the work as possible,
  so that the algorithm should coverge fairly rapidly. The function
  \code{\link{EbuildT}} does the E-step of the EM-algorithm and provides
  more details.

  Note that if there are both missing values in \code{Y} and
  \code{FALSE} entries in \code{Q}, one additional assumption is
  required:  \emph{local independence}.  This assumption states that
  columns of \code{Y} are independent given \code{X}, and hence the
  off-diagonal elements of the covariance matrix is zero.  \code{BQreg}
  does its best to estimate the residual covariance matrix with the
  complete data to test this assumption, but beware of the possible bias
  in the off-diagonal elements, especially with a large fraction of
  missing data.

  
}
\value{
  A list with five components
  \item{B}{A numeric matrix of the same size as \code{Q} giving the regression
    coefficients.}
 \item{b0}{A numeric vector of length \code{nrow{Q}} which gives the
   constant terms.}
 \item{Syy.x}{A symmetric numeric matrix with the same number of rows as
   \code{Q} which provides the covariance of the residuals.}
 \item{converged}{A logical scalar indicating whether or not the EM
   algorithm converged (always \code{TRUE} if \code{Y} is completely
   observed).}
 \item{iterations}{An integer scalar giving the number of EM iterations
   used.  (0 if \code{Y} is completely observed).}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Russell Almond}
\note{

  There are probably way too many bells and whistles on this, but this
  was what I needed for a specific problem.
  
}
\seealso{
  \code{\link{TQB}}, \code{\link{EbuildT}}, \code{\link{matSweep}}
}
\examples{

data(twoskill)  

## X is first two rows, Y is last four.

Est1 <- BQreg(twoskill.data1[,1:2],twoskill.data1[,3:6],
              twoskill.Q)

data(mvMux)

mvX <- mvMux[,1:2] # Complete Data
mvY <- mvMux[,3:5] # Data with missing values.


mvEst <- BQreg(mvX,mvY,trace=TRUE)      # doen't quite converge
mvEst <- BQreg(mvX,mvY,maxstep=50)      # this is enough to converge



}
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
\keyword{ regression }% __ONLY ONE__ keyword per line
